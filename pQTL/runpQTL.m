function runpQTL
clc

%pQTL analysis for rare LoF/LoF+alphamissense masks using the
%burden masks generated by REGENIE under QC = 0 and adjustments of age and
%sex. While adjustment doesn't matter here, it only affects the final
%sample size in the fam files of burden masks --> we pick the largest fam
%file (from a a binary trait with no exclusion criteria).

if isempty(gcp("nocreate"))
    parpool("Processes", 20);
end

qc = 1; % unrelated Europeans (PMID: 39653778)

% read the burden masks ===================================================
bpath = fullfile(fileparts(pwd), "Liver_outcomes", "Results", "Overall");
isoform = getfilenames(bpath, dir=true);

snp = struct;
bed_eid = cell(numel(isoform)*2, 1); % sanity check
for k = 1:numel(isoform)

    % pick the largest sample size
    files = getfilenames(fullfile(bpath, isoform(k)), ".fam", fullpath=true).fam;
    files = arrayfun(@dir, files);
    files = struct2table(files);
    files = sortrows(files, "bytes", "descend");
    
    [~, file_name] = fileparts(string(files.name{1}));
    bed_file = fullfile(bpath, isoform(k), file_name + ".bed");
    bim_file = fullfile(bpath, isoform(k), file_name + ".bim");
    
    % get LoF/LoF+alphamissense masks
    bim_file = bimreader(bim_file, struct=false);

    lof_mask = bim_file.snp(bim_file.snp.endsWith(".mask_lof.0.01"));
    missense_mask = bim_file.snp(bim_file.snp.endsWith("mask_lof_AlphaMissense.0.01"));
    
    snp.chr2.snp(2*k-1, 1) = isoform(k) + "_lof";
    geno = struct;
    [geno.bed, ~, geno.fam] = bedreader(bed_file, lof_mask);
    geno.bed = sum(geno.bed, 2);
    if mean(geno.bed)/2 > 0.5
        geno.bed = 2 - geno.bed;
    end
    [snp.chr2.a1(2*k-1, 1), snp.chr2.a2(2*k-1, 1)] = deal("");
    snp.chr2.chr(2*k-1, 1) = "2";
    snp.chr2.pos(2*k-1, 1) = nan;
    snp.chr2.bed(:, 2*k-1) = geno.bed;
    bed_eid{2*k-1} = geno.fam;

    snp.chr2.snp(2*k, 1) = isoform(k) + "_lofmissense";
    geno = struct;
    [geno.bed, ~, geno.fam] = bedreader(bed_file, missense_mask);
    geno.bed = sum(geno.bed, 2);
    if mean(geno.bed)/2 > 0.5
        geno.bed = 2 - geno.bed;
    end
    [snp.chr2.a1(2*k, 1), snp.chr2.a2(2*k, 1)] = deal("");
    snp.chr2.chr(2*k, 1) = "2";
    snp.chr2.pos(2*k, 1) = nan;
    snp.chr2.bed(:, 2*k) = geno.bed;
    bed_eid{2*k} = geno.fam;

    if k == 1
        snp.chr2.eid = geno.fam;
    end

end

% sanity check on bed eids
bed_eid = horzcat(bed_eid{:});
assert(all(all(bed_eid == bed_eid(:, 1), 1)))


% proteomic files
if isfile("pfiles.mat")
    files = load("pfiles").files;
else

    % find Olink proteins from dictionary file
    dict = load("data_dict_codings_rap.mat").data_dict;
    dict(dict.entity ~= "olink_instance_0", :) = [];
    dict(dict.name == "eid", :) = [];
    dict = dict(1:4, :);
    dict.tag = dict.Field.extractBefore(";");
    phenoParser(query=strings(height(dict), 1), ...
        df=dict.FieldID, ...
        tag=dict.tag, ...
        ukbrap=true, ...
        surv=false, ...
        instance=0, ...
        desc=dict.Field);

    files = matlab.lang.makeValidName(dict.tag);
    save("pfiles.mat", "files")
end

outdir = "raw";
if ~isfolder(outdir), mkdir(outdir), end

% OLINK_BATCH can be defined by running getOLINKbatch function of MAGE.
cvs = ["Age_2", "Age2XSex", "AgeXSex", "Age", "Sex", "BMI", "OLINK_BATCH"];
fid(1:numel(files), 1) = parallel.FevalFuture;
for k = 1:numel(files)

    snps = snp;

    topts = struct;
    topts.trait = files(k);
    topts.output = fullfile(outdir, files(k), files(k));
    if isfile(topts.output + ".xlsx"), continue; end
    topts.qc = qc;
    topts.regenie = false;
    topts.betaround = 9;
    topts.ciround = 6;
    topts.covar = cvs;
    topts.gpu = false;
    topts.catCovar = "OLINK_BATCH";
    topts = namedargs2cell(topts);

    % gwasrunner(snps, topts{:}); %#DEBUG
    fid(k, 1) = parfeval(@gwasrunner, 0, snps, topts{:});

end

if ~all(ismember({fid.State}, 'unavailable'))
    wait(fid)
end

% merge all summary stat files
tab = fileDatastore(fullfile(outdir, files, files + ".xlsx"), ...
    "ReadFcn", @(x)readtable(x, "TextType", "string", ...
    "VariableNamingRule", "preserve"));
tab = gather(tall(tab));
tab = vertcat(tab{:});
tab.P = tab.("P-adj");
tab(:, ["P-adj", "95% CI", "CHR", "POS", "A1", "A2"]) = [];
tab(ismissing(tab.P), :) = [];
tab = renamevars(tab, ["Pheno", "A2FREQ", "β"], ["Trait", "A1Freq", "Beta"]);

idx = tab.A1Freq > 0.5;
tab.A1Freq(idx) = 1 - tab.A1Freq(idx);
tab.Beta(idx) = -tab.Beta(idx);

tab = renamevars(tab, "SNP", "Isoform");
tab = sortrows(tab, "P", "ascend");

% Cauchy P of LoF and LoF/AlphaMissense
for k = 1:numel(isoform)
    df = tab(tab.Scheme.startsWith(isoform(k) + "_"), :);
    usc = unique(df.Scheme);
    df1 = df(df.Scheme == usc(1), :);
    df2 = df(df.Scheme == usc(2), :);

    [f1, f2] = ismember(df1.Trait, df2.Trait); f2(f2<1) = [];
    df1(~f1, :) = [];
    df2 = df2(f2, :);
    p_cauchy = nan(height(df1), 1);
    for j = 1:numel(p_cauchy)
        p_cauchy(j) = CCT([df1.P(j); df2.P(j)]);
    end
    
    p_cauchy_fdr = padjust(p_cauchy, "method", "BH");
    df1.P_cauchy = p_cauchy;
    df2.P_cauchy = p_cauchy;
    df1.P_cauchy_fdr = p_cauchy_fdr;
    df2.P_cauchy_fdr = p_cauchy_fdr;

    df1 = sortrows(df1, "P_cauchy_fdr", "ascend");
    df2 = sortrows(df2, "P_cauchy_fdr", "ascend");

    writetable(df1, "pQTL_ss.xlsx", Sheet=usc(1))
    writetable(df2, "pQTL_ss.xlsx", Sheet=usc(2))

end

end % END
function do_metaanalysis_burden
%@08APR2025: meta-anlysis and ACAT on LoF and LoF+deleterious associations
%for cirrhosis and HCC in UKBB, MVP and Milan studies.

clc

% compile all summary stats
files = struct;
files.name = [fullfile(fileparts(pwd), "Liver_outcomes/Results/rare.xlsx"), ...
    "summary_stat/MVP.xlsx", ...
    "summary_stat/MVP.xlsx",...
    "summary_stat/Milan.xlsx"];
files.study = ["UKBB", "MVP_EUR", "MVP_AFR", "Milan"];
ds = getRareMetaData(files);

% meta-analysis based on traits-mask-isoform
if ~isfile("res.mat")
    res = struct;
    % all studies
    res.three = callMetal(ds);
    
    % repeat without Milan
    res.two = callMetal(ds(ds.Study ~= "Milan", :));
    save("res.mat", "res")
else
    res = load("res.mat").res;
end

% plot and write results --------------------------------------------------
fi = string(fieldnames(res));

for k = 1:numel(fi)
    call_effplotter(res.(fi(k)))
end


end % END

%% subfunctions ===========================================================
function ds = getRareMetaData(files)

for k = 1:numel(files.name)
    if files.study(k) == "UKBB"
        sheet = "Burden";
    elseif files.study(k).contains("MVP")
        sheet = sheetnames(files.name(k));
    else
        sheet = "Sheet1";
    end

    if files.study(k).startsWith("MVP")

        tmp = readtable(files.name(k), Sheet=sheet, TextType="string", VariableNamingRule="preserve");
        if files.study(k).endsWith("_EUR")
            tmp(tmp.Population.lower ~= "european", :) = [];
        else
            tmp(~tmp.Population.lower.contains("african"), :) = [];
        end
        tmp.Population = [];

        tmp = tmp(:, ["Outcome", "Mask", "BETA", "SE", ...
            "P", "N", "Isoform"]);
        tmp = renamevars(tmp, ["BETA", "Outcome"], ["Beta", "Pheno"]);
        tmp{:, ["N_case", "N_control", "MAC_case", "MAC_control"]} = nan;

        idx = tmp.Mask.lower.contains("missense");
        tmp.Mask(idx) = "LoF/Missense";
        tmp.Mask(~idx) = "LoF";
    
    else
        tmp = readtable(files.name(k), Sheet=sheet, TextType="string", VariableNamingRule="preserve");

        if files.study(k) == "UKBB"
            tmp(~ismember(tmp.Pheno, ["HCC", "Cirrhosis"]) | ~ismember(tmp.Mask, ["LoF + AlphaMissense", "LoF"]), :) = [];

            idx = tmp.Mask.contains("AlphaMissense");
            tmp.Mask(idx) = "LoF/Missense";

            tmp.Beta = log(tmp.("Beta|OR(Firth)"));
            tmp = tmp(:, ["Pheno", "Mask", "Beta", "SE", "P", "N", "N.case", "MAC.case", "MAC.control", "Scheme"]);
            tmp = renamevars(tmp, ["MAC.case", "MAC.control", "N.case", "Scheme"], ...
                ["MAC_case", "MAC_control", "N_case", "Isoform"]);
            tmp.N_control = tmp.N - tmp.N_case;
        else % Milan

            tmp(~ismember(tmp.Pheno, ["HCC", "Cirrhosis"]), :) = [];
            tmp.Beta = log(tmp.OR_Firth);
            tmp = tmp(:, ["Pheno", "Group", "n_CASEs", "n_CTRLs", "AC_CASEs", ...
                "AC_CTRLs", "pval_BURDEN", "SE", "Beta"]);

            tmp = renamevars(tmp, ["Group", "n_CASEs", ...
                "n_CTRLs", "AC_CASEs", "AC_CTRLs", "pval_BURDEN"],...
                ["Mask", "N_case", "N_control", "MAC_case", "MAC_control",...
                "P"]);

            idx = tmp.Mask.lower.contains("lof");
            tmp.Isoform = tmp.Mask;
            tmp.Mask(idx) = "LoF";
            tmp.Mask(~idx) = "LoF/Missense";
            tmp.Isoform = tmp.Isoform.erase(["_overall", "- LoF"]).strip;
            
            % Milan's N is variable per isoforms for some reasons: we take max
            gtmp = groupsummary(tmp, "Pheno", @max, ["N_case", "N_control"]);
            gpheno = unique(gtmp.Pheno);
            for g = 1:numel(gpheno)
                fidx = tmp.Pheno == gpheno(g);
                tmp.N_case(fidx) = gtmp.fun1_N_case(gtmp.Pheno == gpheno(g));
                tmp.N_control(fidx) = gtmp.fun1_N_control(gtmp.Pheno == gpheno(g));
            end

            tmp.N = tmp.N_control + tmp.N_case;

        end
    end

    tmp.Study(:) = files.study(k);
    files.res{k, 1} = tmp;
end

ds = vertcat(files.res{:});
ds.Isoform = ds.Isoform.upper;

idx = ds.Isoform == "APOB48100";
ds.Isoform(idx) = "APOB48/100";

ds.Isoform = ds.Isoform.replace("APOB", "ApoB");

end % END

%% ========================================================================
function res = callMetal(df)

df.id = df.Pheno + ":" + df.Mask + ":" + df.Isoform;

% sanity check
acheck = groupsummary(df, "Study", @(x)numel(duplicates(x)), "id");
assert(all(acheck.fun1_id == 0))

us = unique(df.Study);
for j = 1:numel(us)
    tmp = df(df.Study == us(j), :);

    % dummy allele/freq
    tmp.ALLELE1(:) = "A";
    tmp.ALLELE0(:) = "T";
    tmp.A1FREQ(:) = 0.3;
    tmp.CHROM(:) = 2;
    tmp.GENPOS(:) = 200;
    tmp(isinf(tmp.SE) | ismissing(tmp.SE), :) = [];
    writetable(tmp, us(j) + ".txt", "Delimiter", ",");
    clear tmp
    dos2unix(fullfile(pwd, us(j) + ".txt"), verbose=true);
end

metal(us + ".txt", snp="id", ea="ALLELE1", nea="ALLELE0", eaf="A1FREQ", ...
    effect="Beta", n="N", p="P", se="SE", GENOMICCONTROL="OFF", ...
    OUTFILE="meta", chr="CHROM", pos="GENPOS", SCHEME="STDERR", ...
    SEPARATOR="COMMA")

tmp = readtable("meta.txt", ...
    TextType="string", ...
    VariableNamingRule="preserve");
tmp = renamevars(tmp, "MarkerName", "id");

% add pheno/isoform/mask
meta_data = tmp.id.split(":");
tmp.Pheno = meta_data(:, 1);
tmp.Mask = meta_data(:, 2);
tmp.Isoform = meta_data(:, 3);
tmp = tmp(:, ["Isoform", "Pheno", "Mask", "Effect", ...
    "StdErr", "P-value","N", "Direction", "HetISq", "HetChiSq", ...
    "HetDf", "HetPVal", "id"]);
tmp = renamevars(tmp, ["Effect", "StdErr", "P-value"], ["Beta", "SE", "P"]);
tmp.Study(:) = us.join(",");

% add sum N_case/control and MAC_case/control
ncols = ["N_case", "N_control", "MAC_case", "MAC_control"];
nsum = groupsummary(df, "id", @sum, ncols);
nsum.Properties.VariableNames = erase(colnames(nsum), "fun1_");
[~, ff] = ismember(tmp.id, nsum.id);
for n = 1:numel(ncols)
    tmp.(ncols(n)) = nsum.(ncols(n))(ff);
end
tmp.id = [];

% add each study's summary stat to pooled results
df.Direction(:) = "";
df{:, ["HetChiSq", "HetDf", "HetISq", "HetPVal"]} = nan;

res = [tmp; df];

clear tmp df

del_files = ["meta.txt.info", "meta.txt", us' + ".txt"];
arrayfun(@delete, del_files)


% ACAT method for each trait-isoform pair  --------------------------------
res.grp = res.Pheno + ":" + res.Isoform + ":" + res.Scheme;

res.ACAT_P(:) = nan;

ugr = unique(res.grp);

for j = 1:numel(ugr)
    idx = res.grp == ugr(j);
    tab = res(idx, :);
    tab(~tab.Study.contains(","), :) = []; % remove single studies

    assert(height(tab) == numel(unique(tab.Mask)))

    % p-value combination method using the Cauchy distribution
    res.ACAT_P(idx) = CCT(tab.P);
end

res.grp = [];

end % END

%% ========================================================================
function call_effplotter(ds)
% plot meta-analysis results.

% separate output directories for 2 and 3 studies 
out_dir = unique(ds.Study);
out_dir(~out_dir.contains(",")) = [];
out_dir = out_dir.replace(",", "_");
out_dir = fullfile(pwd, out_dir);

if ~isfolder(out_dir), mkdir(out_dir); end 

% plots per each strata: mask + pheno
ds.str = ds.Pheno + ":" + ds.Mask;

writetable(ds, fullfile(out_dir, "meta_res.xlsx"))
ustr = unique(ds.str);

for k = 1:numel(ustr)
    
    out_name = fullfile(out_dir, matlab.lang.makeValidName(ustr(k)));
    if isfile(out_name + ".pdf"), continue; end
    df = ds(ds.str == ustr(k), :);

    df.Study = df.Study.replace("_", "-");
    
    fe_idx = df.Study.contains(",");
    df.Study(fe_idx) = "Fixed-effect";
    df.Marker(:) = "square";
    df.Marker(fe_idx) = "diamond";
    df.N(df.Study == "Fixed-effect") = max(df.N(df.Study == "Fixed-effect")); % maximum sample size

    df.Study = df.Study + ", N = " + df.N;
    
    df.OR = exp(df.Beta); df.Beta = [];
    
    marker_tab = struct;
    marker_tab.a1 = unique(df.Study);
    marker_tab.a2 = repmat("O", numel(marker_tab.a1), 1);
    marker_tab.a2(marker_tab.a1.lower.contains("fixed-effect")) = "diamond";
    marker_tab = struct2table(marker_tab);

    idx = isinf(df.SE) | ismissing(df.SE);
    df{idx, ["OR", "SE", "P"]} = nan;

    effPlotter(df, ...
        betaColumn="OR", ...
        binary=true(height(df), 1), ...
        groupCol="Study", ...
        yCol="Isoform", ...
        save=true, ...
        hide=true, ...
        squeeze=true, ...
        yOrder=["ApoB", "ApoB100", "ApoB48/100"], ...
        marker=marker_tab, ...
        sortEffect=false, ...
        markersize=200, ...
        fontsize=21, ...
        ciLineWidth=2, ...
        format=["pdf", "png"], ...
        output=out_name, ...
        ybold=true, ...
        xlabel="OR", ...
        xbold=true, ...
        xscale="log", ...
        squeezeOffset=eps, ...
        legFontOffset=4, ...
        colormap="inferno")

end

end % END
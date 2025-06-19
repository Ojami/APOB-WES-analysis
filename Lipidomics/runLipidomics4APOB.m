function runLipidomics4APOB
%@17OCT2024: rare variant association with lipidomics data.

pheno = readtable("NMR_Metabolomics_meta.xlsx", ...
    TextType="string", ...
    VariableNamingRule="preserve");

% get path to pheno/covar/pred files from REGENIE step 1 ------------------
phewide_path = fullfile(pwd, "step1");
pheno.traits = matlab.lang.makeValidName(pheno.name_ukb);

for k = 1:numel(pheno.traits)
    pheno.pred(k) = getfilenames(fullfile(phewide_path, pheno.traits(k)), ...
        "list", fullpath=true).list;
    files = getfilenames(fullfile(phewide_path, pheno.traits(k)), ...
        "txt", fullpath=true).txt;
    pheno.covarFile(k) = files(files.endsWith(".covarFile.txt"));
    pheno.phenoFile(k) = files(files.endsWith(".phenoFile.txt"));
end

varset = load(fullfile(fileparts(pwd), "APOB.varset.mat"));
vep = varset.vep;
lc_lof_ids = vep.id(vep.LOFTEE == "LC"); % to remove based on LOFTEE predictions
varset = varset.apob;
fi = string(fieldnames(varset));

threads = 10;
genes = "APOB";
aaf_bins = 0.01; 

for k = 1:numel(fi) % loop over isoforms

    rarevariants = setdiff(varset.(fi(k)), lc_lof_ids);

    if ~isfolder(fi(k))
        gwasCaller(genes=genes, ...
            report=false, ...
            workers=threads, ...
            pheno_genes=repmat({pheno.traits},1,numel(genes)), ...
            pred=repmat({pheno.pred}, 1, numel(genes)), ...
            phenoFile=repmat({pheno.phenoFile}, 1, numel(genes)), ...
            covarFile=repmat({pheno.covarFile}, 1, numel(genes)),...
            ignore_gene_scan=true, ...
            ignore_regional_plots=true, ...
            gene_scan_phewas=false, ...
            dir=fi(k), ...
            varset=rarevariants, ...
            aaf_bins=aaf_bins, ...
            weshome="/dnax/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/")
    end

end

sheet = "0.01"; % MAF < 0.01
test = "ADD"; % Burden test
pth = fullfile(fileparts(pwd), "CircularPlot"); % R code is there
tmp = readtable(fullfile(pth, "template.xlsx"), ...
    TextType="string", ...
    VariableNamingRule="preserve");

% create input tables for circle plots (Carolin's pipeline)
for k = 1:numel(fi)
    out = tmp;

    file = fullfile(pwd, fi(k), "APOB", "APOB.RVA.xlsx");
    df = readtable(file, TextType="string", VariableNamingRule="preserve",...
        Sheet=sheet);
    df = df(:, ["Pheno", "Mask", "β|OR(Firth)", "SE", test]);
    
    if fi(k) == "RVAs"
        mask_name = "mask_lof_REVEL5ORCAD20.0.01";
    else
        mask_name = "mask_lof.0.01";
    end

    df(df.Mask ~= mask_name, :) = [];
    
    % add Estimate, StdErr, p.value
    [f1, f2] = ismember(out.name_ukb, df.Pheno); 
    out(~f1, :) = []; df = df(f2(f1), :);
    out.Estimate = df.("β|OR(Firth)");
    out.StdErr = df.SE;
    out.("p.value") = df.(test);
    

    writetable(out, fullfile(pth, fi(k) + ".xlsx"))
    clear out

end

end % END
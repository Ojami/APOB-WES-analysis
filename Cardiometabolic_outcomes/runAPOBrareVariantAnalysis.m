function runAPOBrareVariantAnalysis
clc
% performs gene-based rare variant analysis with REGENIE for APOB rare
% variants with the following adjustments (Fig. 2):
% age, sex, age×sex, age2×sex, age2, smoking, alcohol consumption, with
% additional adjustments for diabetes and hypertension for coronary artery
% disease.

%% run REGENIE step 1 to fit models ---------------------------------------
traits = struct;
traits.base = ["ALT", "LDL", "HDLCholesterol", "Triglycerides", "Glucose", ...
    "HbA1c", "C_reactiveProtein", "BMI", "ChronicKidneyFailure", ...
    "HeartFailure_HF_", "CoronaryArteryDisease_CAD_"];
traits.mri = "PDFF";

cvr = struct;
cvr.base_all = ["Age", "Sex", "BMI", "SmokingStatus", "AlcoholConsumption_g_day_"];
cvr.base_cat = "SmokingStatus";
cvr.mri_all = ["AgeAtMRI", "Sex", "BMI_MRI", "SmokingStatusMRI", "AlcoholConsumption_g_day_MRI"];
cvr.mri_cat = "SmokingStatusMRI";

step1_path = fullfile(pwd, "step1");
if ~isfolder(step1_path), mkdir(step1_path); end

for k = 1:numel(traits.base)

    out_dir = fullfile(step1_path, traits.base(k));
    if ~isfolder(out_dir), mkdir(out_dir); end
    
    % remove covariate if it's the response variable
    tmp_cvr = cvr.base_all;
    tmp_cvr = setdiff(tmp_cvr, traits.base(k));

    if traits.base(k) == "CoronaryArteryDisease_CAD_"
        tmp_cvr = union(tmp_cvr, ["Diabetes", "Hypertension"]);
    end

    call_regenie(step=1, ...
        bed="/dnax/GRM/ukb.qcn.indep.chr1-22.bed", ...
        trait=traits(k), ...
        covar=tmp_cvr, ...
        agesexInt=true, ...
        threads=30, ...
        bsize=1000,...
        out=fullfile(out_dir, traits.base(k)), ...
        catCovarList=cvr.base_cat, ...
        firth=true, approx=true, firth_se=true)
end

for k = 1:numel(traits.mri)

    out_dir = fullfile(step1_path, traits.mri(k));
    if ~isfolder(out_dir), mkdir(out_dir); end

    % remove covariate if it's the response variable
    tmp_cvr = cvr.mri_all;
    tmp_cvr = setdiff(tmp_cvr, traits.mri(k));

   call_regenie(step=1, ...
        bed="/dnax/GRM/ukb.qcn.indep.chr1-22.bed", ...
        trait=traits.mri(k), ...
        covar=tmp_cvr, ...
        agesexInt=true, ...
        threads=30, ...
        bsize=1000,...
        out=fullfile(out_dir, traits.mri(k)), ...
        catCovarList=cvr.base_cat, ...
        firth=true, approx=true, firth_se=true)
end

%% run rare variant analysis for each subset of APOB variants -------------
aaf_bins = 0.01; 
threads = 10;
genes = "APOB";

varset = load(fullfile(fileparts(pwd), "APOB.varset.mat"));
vep = varset.vep;
lc_lof_ids = vep.id(vep.LOFTEE == "LC"); % to remove based on LOFTEE predictions
varset = varset.apob;
fi = string(fieldnames(varset)); % APOB isoforms 

traits = union(traits.base, traits.mri);

[covarFile, phenoFile, predfiles, catCovar] = deal(strings(numel(traits), 1));
for k = 1:numel(traits)
    tmp = getfilenames(fullfile(step1_path, traits(k)), ".list", fullpath=true).list;
    predfiles(k) = tmp(tmp.endsWith("_pred.list"));

    tmp = getfilenames(fullfile(step1_path, traits(k)), ".txt", fullpath=true).txt;
    covarFile(k) = tmp(tmp.endsWith(".covarFile.txt"));
    phenoFile(k) = tmp(tmp.endsWith(".phenoFile.txt"));

    cat_covar = bfilereader(covarFile(k), "summary", "firstline");
    catCovar(k) = cat_covar(cat_covar.lower.contains("smoking"));

end


for k = 1:numel(fi) % loop over isoforms (REGENIE step 2)
    rarevariants = setdiff(varset.(fi(k)), lc_lof_ids);
    
    pth_out = fullfile(pwd, "Results", "Overall", fi(k));
    if isfolder(pth_out), continue; end
    
    % gwasCaller by default activates firth correction before sending the
    % jobs to REGENIE
    gwasCaller(genes=genes, ...
        report=false, ...
        workers=threads, ...
        pheno_genes=repmat({traits},1,numel(genes)), ...
        pred=repmat({predfiles}, 1, numel(genes)), ...
        covarFile=repmat({covarFile}, 1, numel(genes)), ...
        phenoFile=repmat({phenoFile}, 1, numel(genes)), ...
        catCovarCell=repmat({catCovar}, 1, numel(genes)), ...
        ignore_gene_scan=true, ...
        ignore_regional_plots=true, ...
        gene_scan_phewas=false, ...
        varset=rarevariants, ...
        aaf_bins=aaf_bins, ...
        dir=pth_out, ...
        weshome="/dnax/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/")

end


end % END
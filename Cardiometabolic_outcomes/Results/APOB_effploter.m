function APOB_effploter
% plots the summary statistics from REGENIE step 2 to tables
% run this helper function after running runAPOBrareVariantAnalysis
% function.
clc

traits = "";

out_dir = fullfile(pwd, "plots");
if ~isfolder(out_dir), mkdir(out_dir); end

pth = fullfile(pwd, "Overall");
di = getfilenames(pth, dir=true);
files = struct;
files.di = di;
files.sc = di;
files = struct2table(files);
files.path = fullfile(pth, files.di, "APOB");

% read all files
ds = struct;
for k = 1:height(files)
    [ds.rare{k, 1}, ds.joint{k, 1}] = fetchTable(files(k, :), traits);
end
ds.rare(cellfun(@isempty, ds.rare)) = [];
ds.joint(cellfun(@isempty, ds.joint)) = [];
ds.rare = vertcat(ds.rare{:});
ds.joint = vertcat(ds.joint{:});

% effect plots
fis = string(fieldnames(ds));
for k = 1:numel(fis)
    
    tab = ds.(fis(k));
    if isempty(tab), continue; end
    if any(tab.test == "GENE-P")
        joint_flag = true;
        tab(tab.test ~= "GENE-P", :) = [];
    else
        joint_flag = false;
        tab(tab.test ~= "Burden", :) = [];

        % @18APR2025: add Gene-P to the secondary y-axis
        jtab = ds.joint;
        jtab(jtab.test ~= "GENE-P", :) = [];
        jtab.tmp = jtab.Pheno + ":" + jtab.bin;
        tab.tmp = tab.Pheno + ":" + tab.bin;
        tab.y2(:) = "";
        [~, ff] = ismember(tab.tmp, jtab.tmp);
        tab.y2 = compose("%.2E", jtab.P(ff));
        clear jtab
    end


    fig_name = fullfile(out_dir, fis(k));

    if joint_flag
        tab.beta = -log10(tab.P);
        tab(ismember(tab.Pheno, ["Triglycerides", "LDL", "Cholesterol"]), :) = [];
        
        effPlotter(tab, ...
            betaColumn="beta", ...
            groupCol="bin", ...
            ignoreCI=true,...
            save=true, ...
            yCol="Pheno", ...
            hide=true, ...
            squeeze=true, ...
            xlabel="-log_{10} GENE-P", ...
            visible="off", ...
            legTitle="Isoform", ...
            colormap="magma", ...
            fontsize=13, ...
            sortEffect=false, ...
            output=fig_name, ...
            binary=false(height(tab), 1))

    else
        % convert OR to beta
        or_idx = ~ismissing(tab.("A1F.case"));
        tab.beta(or_idx) = log(tab.beta(or_idx));

        tab(ismember(tab.Pheno, ["Cholesterol", "Cardiovascular disease (CVD)", "Hypertension"]), :) = [];
        tab.Pheno(tab.Pheno == "LDL") = "LDL-C";
        tab.Pheno(tab.Pheno == "HDL cholesterol") = "HDL-C";
        tab.Pheno(tab.Pheno == "Glycated haemoglobin (HbA1c)") = "HbA1c";
        tab.Pheno(tab.Pheno.contains("(HF)")) = "Heart failure";
        tab.Pheno(tab.Pheno.contains("(CAD)")) = "Coronary artery disease";
        
        yOrder = ["LDL-C", "HDL-C", "Triglycerides", ...
            "Coronary artery disease", "Heart failure",...
            "ALT", "PDFF",...
            "Chronic kidney failure", "Glucose", "HbA1c",...
            "Diabetes", "C-reactive protein", "BMI"]; 

        % with Gene-P printed   
        % effPlotter(tab, betaColumn="beta", groupCol="bin", ignoreCI=false,...
        %     save=true, yCol="Pheno", ...
        %     hide=true, squeeze=true, ...
        %     xlabel="Beta (burden)", visible="off", ...
        %     legTitle="Isoform", colormap="colorblind", ...
        %     fontsize=13, sortEffect=false, ...
        %     output=fig_name, ...
        %     binary=false(height(tab), 1), yOrder=yOrder, ...
        %     format=["png", "pdf"], ycol2="y2", ycol2Label="Gene-P", ...
        %     fontsizeYcol2=10)
        
        effPlotter(tab, ...
            betaColumn="beta", ...
            groupCol="bin", ...
            ignoreCI=false,...
            save=true, ...
            yCol="Pheno", ...
            hide=true, ...
            squeeze=true, ...
            xlabel="Beta (burden)", ...
            visible="off", ...
            legTitle="Isoform", ...
            colormap="colorblind", ...
            fontsize=16, ...
            sortEffect=false, ...
            output=fig_name, ...
            binary=false(height(tab), 1), ...
            yOrder=yOrder, ...
            format=["png", "pdf"], ...
            fontsizeYcol2=10, ...
            markersize=90, ...
            ybold=true, ...
            xbold=true)
        
    end
end

end % END

%% subfunctions ===========================================================
function [rtab, jtab] = fetchTable(di, traits)

file = getfilenames(di.path, "xlsx", fullpath=true).xlsx;

rtab = readtable(file, TextType="string", VariableNamingRule="preserve", Sheet="0.01");
jtab = readtable(file, TextType="string", VariableNamingRule="preserve", Sheet="Joint");

if all(traits == "")
    traits = unique(rtab.Pheno);
end

% keep only LoF mask and custom traits
rtab(~ismember(rtab.Pheno, traits) | rtab.Mask ~= "mask_lof.0.01", :) = []; 
rtab.Mask(:) = "LoF";

jtab(~ismember(jtab.Pheno, traits), :) = [];
rtab = renamevars(rtab, "β|OR(Firth)", "beta"); % beta or OR (Firth)

spheno = setdiff(traits, rtab.Pheno);
if ~isempty(spheno)
    tmp = repmat(rtab(1, :), numel(spheno), 1);
    tmp.Pheno = spheno;
    tmp{:, ["beta", "SE", "ADD", "ADD-ACATO", "ADD-ACATV",...
        "ADD-SKAT", "ADD-SKATO", "ADD-SKATO-ACAT"]} = nan;
    rtab = [rtab; tmp];
end

% modify col names
cols = colnames(rtab) == "ADD";

rtab = renamevars(rtab, cols, "Burden");
rtab.Properties.VariableNames = erase(colnames(rtab), textBoundary("start") + "ADD-");
jtab.Properties.VariableNames = erase(colnames(jtab), textBoundary("start") + "ADD-");
jtab = renamevars(jtab, ...
    ["GENE_P", "BURDEN-SBAT_NEG", "BURDEN-SBAT_POS"], ...
    ["GENE-P", "BURDEN-SBAT_{Neg}", "BURDEN-SBAT_{Pos}"]);

% stack
rtab = stack(rtab, ...
    ["Burden", "ACATO", "SKAT", "ACATV", "SKATO", "SKATO-ACAT"], ...
    IndexVariableName="test", NewDataVariableName="P");
jtab = stack(jtab,...
    ["ACATV-ACAT", "BURDEN-ACAT", "BURDEN-SBAT", "BURDEN-SBAT_{Neg}", ...
    "BURDEN-SBAT_{Pos}", "GENE-P", "SKATO-ACAT"], ...
    IndexVariableName="test", NewDataVariableName="P");
rtab.test = string(rtab.test);
jtab.test = string(jtab.test);

rtab.P(rtab.P == 0) = realmin;
jtab.P(jtab.P == 0) = realmin;

% add Beta
beta_col = unique(rtab.test); 
bet_col = beta_col(floor(numel(beta_col)/2) + 1);
rtab.y2 = compose("%.3G", rtab.beta);
rtab.y2(rtab.test ~= bet_col) = missing;

if any(colnames(di) == "q")
    rtab.bin(:) = di.q;
    jtab.bin(:) = di.q;
    di.sc = di.sc.replace("RVA", "ApoB");
    di.sc = di.sc.replace(["a48", "b48"], ["48/100", "100"]);
    [jtab.Scheme(:), rtab.Scheme(:)] = deal(di.sc);
else
    di.di = di.di.replace("RVA", "ApoB");
    di.di = di.di.replace(["a48", "b48"], ["48/100", "100"]);
    rtab.bin(:) = di.di;
    jtab.bin(:) = di.di;
end

end % END

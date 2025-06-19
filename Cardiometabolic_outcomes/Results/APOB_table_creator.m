function APOB_table_creator
% writes summary statistics from REGENIE step 2 to tables
% run this helper function after running runAPOBrareVariantAnalysis
% function.
clc

traits = ""; % fetch all traits

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
    [ds.rare{k, 1}, ds.joint{k, 1}, ds.singleton{k, 1}] = fetchTable(files(k, :), traits);
end
ds.rare(cellfun(@isempty, ds.rare)) = [];
ds.joint(cellfun(@isempty, ds.joint)) = [];
ds.singleton(cellfun(@isempty, ds.singleton)) = [];
ds.rare = vertcat(ds.rare{:});
ds.joint = vertcat(ds.joint{:});
ds.singleton = vertcat(ds.singleton{:});

% saving and plotting
fis = string(fieldnames(ds));
for k = 1:numel(fis)

    tab = ds.(fis(k));
    utest = unique(tab.test);

    if any(colnames(tab) == "beta")
        tab = renamevars(tab, "beta", "Beta|OR(Firth)");
    end
    tab = renamevars(tab, "bin", "Scheme");   
    
    % write to table 
    if ~isfile(fis(k) + ".xlsx")
        for j = 1:numel(utest)
            tmp = tab(tab.test == utest(j), :);
            tmp(:, ["ID", "CHR", "test"]) = [];
            if any(colnames(tmp) == "y2"), tmp.y2 = []; end
            
            writetable(tmp, fis(k) + ".xlsx", Sheet=utest(j))
        end
    end
    
end

end % END

%% subfunctions ===========================================================
function [rtab, jtab, stab] = fetchTable(di, traits)

file = getfilenames(di.path, "xlsx", fullpath=true).xlsx;

rtab = readtable(file, TextType="string", VariableNamingRule="preserve", Sheet="0.01");
stab = readtable(file, TextType="string", VariableNamingRule="preserve", Sheet="singleton");
jtab = readtable(file, TextType="string", VariableNamingRule="preserve", Sheet="Joint");

if all(traits == "")
    traits = unique(rtab.Pheno);
end

% keep only LoF mask and custom traits
rtab(~ismember(rtab.Pheno, traits) | ~ismember(rtab.Mask, ["mask_lof.0.01" "mask_lof_AlphaMissense.0.01"]), :) = []; 

idx = rtab.Mask.lower.contains("alphamissense");
rtab.Mask(idx) = "LoF + AlphaMissense";
rtab.Mask(~idx) = "LoF";

stab(~ismember(stab.Pheno, traits) | ~ismember(stab.Mask, ["mask_lof.singleton", "mask_lof_AlphaMissense.singleton"]), :) = []; 
idx = stab.Mask.lower.contains("alphamissense");
stab.Mask(idx) = "LoF + AlphaMissense singleton";
stab.Mask(~idx) = "LoF singleton";

jtab(~ismember(jtab.Pheno, traits), :) = [];
rtab = renamevars(rtab, "β|OR(Firth)", "beta"); % beta or OR (Firth)
stab = renamevars(stab, "β|OR(Firth)", "beta"); % beta or OR (Firth)

spheno = setdiff(traits, rtab.Pheno);
if ~isempty(spheno)
    tmp = repmat(rtab(1, :), numel(spheno), 1);
    tmp.Pheno = spheno;
    tmp{:, ["beta", "SE", "ADD", "ADD-ACATO", "ADD-ACATV", ...
        "ADD-SKAT", "ADD-SKATO", "ADD-SKATO-ACAT"]} = nan;
    rtab = [rtab; tmp];
end

spheno = setdiff(traits, stab.Pheno);
if ~isempty(spheno)
    tmp = repmat(stab(1, :), numel(spheno), 1);
    tmp.Pheno = spheno;
    tmp{:, ["beta", "SE", "ADD"]} = nan;
    stab = [stab; tmp];
end


% modify col names
cols = colnames(rtab) == "ADD";
rtab = renamevars(rtab, cols, "Burden");
cols = colnames(stab) == "ADD";
stab = renamevars(stab, cols, "Burden");

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
stab = stack(stab, "Burden", ...
    IndexVariableName="test", NewDataVariableName="P");

rtab.test = string(rtab.test);
jtab.test = string(jtab.test);
stab.test = string(stab.test);

% add Beta
beta_col = unique(rtab.test); 
bet_col = beta_col(floor(numel(beta_col)/2) + 1);
rtab.y2 = compose("%.3G", rtab.beta);
rtab.y2(rtab.test ~= bet_col) = missing;

if any(colnames(di) == "q")
    rtab.bin(:) = di.q;
    jtab.bin(:) = di.q;
    stab.bin(:) = di.q;

    di.sc = di.sc.replace("RVA", "ApoB");
    di.sc = di.sc.replace(["a48", "b48"], ["48/100", "100"]);

    [jtab.Scheme(:), rtab.Scheme(:), stab.Scheme(:)] = deal(di.sc);
else
    di.di = di.di.replace("RVA", "ApoB");
    di.di = di.di.replace(["a48", "b48"], ["48/100", "100"]);

    rtab.bin(:) = di.di;
    jtab.bin(:) = di.di;
    stab.bin(:) = di.di;
end

end % END

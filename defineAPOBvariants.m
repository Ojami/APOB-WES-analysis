function defineAPOBvariants
% using the pre-annotated 500 K WES UKBB data using VEP and REGENIE
% anno/set/mask files (see runVEP and geneWrapper functions of MAGE to
% generated these files; alternatively use DNAnexus UKB-RAP preannotated
% data already available on the project), we extract the variants on APOB
% gene and define isoforms.
% Note that APOB is on negative strand: larger aa <--> smaller bp

% read set file for chr 2
anno_path = "/dnax/geneWrapper/annotation.files.500K";
vset = readtable(fullfile(anno_path, "geneWrapper.set.chr2.txt"),...
    TextType="string", ...
    NumHeaderLines=0);

vset(~vset.(1).startsWith("APOB("), :) = [];
vset = split(vset.(4), ",");

apob48End = 21010408; % GRCh38

% a temp bim struct for VEP API to get LOFTEE predictions
bim = struct;
bim.snp = vset;
ids = bim.snp.split(":");
bim.chr = ids(:, 1);
bim.a1 = ids(:, 4);
bim.a2 = ids(:, 3);
bim.pos = double(ids(:, 2));

% define APOB isoforms =====================================================
apob = struct;
apob.RVA = vset; % full APOB

% before/after apob48: APOB is on negative strand
idx = bim.pos < apob48End;
apob.RVAb48 = bim.snp(idx); % 100-only
apob.RVAa48 = bim.snp(~idx); % 48/100

% call VEP API
vep = callVEP(bim, backup=false,...
    verbose=true, ...
    genomeRef="38", ...
    loftee=true, ...
    method="vep");

save("APOB.varset.mat", "apob", "vep")

end % END
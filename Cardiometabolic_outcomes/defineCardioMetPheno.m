function defineCardioMetPheno
% defining cardiometabolic and liver traits (Fig. 2) for rare variant
% analysis with following covariates:
% age, sex, age×sex, age2×sex, age2, smoking, alcohol consumption, with
% additional adjustments for diabetes and hypertension for coronary artery
% diseas
% Notes:
%   1-to define alcohol consumption, see getUKB_alcohol.m function of MAGE
%     toolbox.
%   2-to define age at MRI visit (for PDFF covariate), see
%   getAgeAtInstance.m function of MAGE toolbox.
%   3-see PMID 39653778 for definition of other cardiometabolic traits.

% get continuous traits
qt = struct;
qt.df = ["30620", "30780", "30760", "30870", "30740", "30750", "30710", ...
    "21001", "40061", "21001"];
qt.tag = ["ALT", "LDL", "HDL cholesterol", "Triglycerides", "Glucose", ...
    "HbA1c", "C-reactive protein", "BMI", "PDFF", "BMI_MRI"];
n = numel(qt.df);
qt.instance = [repmat("0", 1, n-2), "2", "2"]; % PDFF: MRI-visit
a = cell(n, 1);
for k = 1:n
    a{k} = phenoParser(query="", ...
        df=qt.df(k), ...
        tag=qt.tag(k), ...
        instance=qt.instance(k), ...
        surv=false, ...
        all=true);
end

% define binary traits (PMID 39653778)
defDiabetes

phenoParser(query=["1074;1075", "I20;I200;I201;I208;I209;I21;I210;I211;I212;I213;I214;I219;I21X;I22;I220;I221;I228;I229;I23;I230;I231;I232;I233;I234;I235;I236;I238;I24;I240;I241;I248;I249;I25;I250;I251;I252;I253;I254;I255;I256;I258;I259"],...
    surv=false, ...
    df=["20002", ""],...
    coding=["", "19"],...
    tag="Coronary artery disease (CAD)", ...
    merge=true, ...
    all=true);

phenoParser(query=["1076", "I50;I500;I501;I509"],...
    surv=false, ...
    df=["20002", ""], ...
    coding=["", "19"], ...
    tag="Heart failure (HF)", ...
    merge=true, ...
    all=true);

phenoParser(query=["1192;1193", "N18;N180;N181;N182;N183;N184;N185;N188;N189"],...
    surv=false, ...
    df=["20002", ""], ...
    coding=["", "19"], ...
    tag="Chronic kidney failure", ...
    merge=true, ...
    all=true);

phenoParser(query=["1065;1072", "I10", "2", ...
    "1140860736;1140866724;1140860422;1140860752;1140875868;1140856568;1140881716;1140861090;1140881702;1140923572;1140861138;1140926780;1140866396;1140866244;1140910620;1140866354;1140888512;1140866422;1140866426;1140866222;1140909722;1140879802;1141200400;1140866704;1140861110;1141145658;1141150898;1140861136;1140917428;1140888578;1140864410;1140866128;1141153006;1140866764;1140866212;1140864550;1140866738;1141146126;1141194810;1141180778;1141146124;1141146128;1140860426;1140881722;1140851556;1140866122;1140866450;1141194794;1141194800;1140866226;1140866546;1140866784;1140866132;1140860356;1140866692;1140916342;1140860304;1140866782;1140851492;1140860266;1140879758;1140860380;1141168964;1141175224;1141184324;1140879760;1140864950;1140860382;1140851360;1140861130;1141187094;1141156836;1140910630;1140860750;1140860764;1140864910;1140860706;1141190934;1140861176;1141171152;1140927934;1140866712;1141199858;1140879822;1140863724;1140909368;1140879762;1140851332;1140866440;1140866138;1140909706;1140866144;1140864202;1140860882;1140923402;1140923276;1140860386;1141201040;1140923404;1140923336;1140923272;1141181186;1141172686;1140866554;1140860194;1140860390;1140860802;1141180598;1141179974;1140916362;1141151018;1140866156;1140875798;1140875870;1140851418;1140875796;1141157136;1140879806;1140926778;1140861166;1141145668;1140851400;1140881894;1140866110;1141188636;1140922714;1141169516;1140866182;1140866402;1140866404;1141150328;1140860492;1140888552;1140860790;1140851338;1141201244;1141171336;1140866164;1140851362;1141168498;1140888646;1141165470;1141145870;1140888556;1140866192;1141167108;1140851414;1140866406;1140928624;1140866116;1140866506;1140866408;1140851412;1140866194;1141169088;1140909708;1141195258;1140866484;1141152600;1140866802;1141152076;1141156754;1140866460;1140866798;1140866800;1141180238;1140866074;1140866162;1140866072;1140866168;1140866146;1140851364;1141151382;1140851432;1141188936;1141167758;1141164148;1140866078;1140866804;1140860394;1140860396;1140860776;1140860784;1141201250;1141152998;1141172682;1140861190;1140861202;1140866410;1140860398;1141150560;1141188920;1141187962;1140860232;1141181520;1140879824;1140860244;1140861276;1140851420;1140866412;1140860400;1140866248;1140881728;1140866334;1141153026;1140879826;1140851784;1140860696;1140864952;1141199940;1140917076;1140860274;1140860402;1140916356;1141151016;1140866084;1140860278;1140917452;1140866094;1140875808;1140866092;1140879818;1140860308;1140860404;1140851522;1141172492;1141187790;1140866220;1140860406;1140866416;1140866420;1140923712;1140860434;1140864176;1140860192;1140860312;1140860316;1141194804;1140916870;1141146378;1140866352;1141164280;1141164276;1140866136;1140866446;1140879810;1140861088;1141157140;1141150538;1140911088;1140861114;1141169730;1140872568;1140926966;1140872472;1140888922;1141162546;1141193282;1141193346;1140917068;1140879830;1140851484;1140868036;1141201814;1140879834;1140860320;1140888560;1141180592;1140860292;1140860322;1140866210;1140928212;1140866102;1140866230;1140861194;1140860410;1140910614;1140866766;1141156808;1140879842;1140860418;1140860738;1140860728;1140860806;1141200698;1141187048;1140851658;1140881712;1140866262;1140866726;1140866466;1141150500;1140916730;1140911698;1140860362;1140879854;1140860332;1140860318;1140866308;1140864574;1140866312;1140866232;1140866318;1140866236;1140851508;1140866306;1140860878;1140851430;1140928234;1141164154;1141153316;1141166006;1141187788;1141146184;1140860352;1140860358;1140860324;1140860328;1140866756;1140927940;1141182968;1141167822;1141171344;1140861128;1140879866;1140860340;1140860342;1141194808;1140860336;1140860330;1140888496;1140864874;1140860172;1140916628;1140860250;1140860904;1141153328;1140860222;1140860334;1140866328;1140866360;1140866388;1140866324;1141180772;1140866330;1140866332;1141195254;1141165476;1141188408;1140881692;1140851336;1140857140;1141145660;1141201038;1140851790;1140866758;1141190160;1140851436;1141187774;1140888510;1141150926;1141151474;1140860338;1140860294;1141187780;1140866108;1141153032;1141174684;1141167832;1140864618;1140860714;1140922324;1141171804"], ...
    surv=false,...
    df=["20002", "", "6153;6177", "20003"], ...
    coding=["", "19", "", ""], ...
    tag="Hypertension",...
    merge=true,...
    all=true);

% Somking status: never/previous/current. Prefer not to answer --> missing
tmp = phenoParser(query="", ...
    df="20116", ...
    surv=false, ...
    all=true, ...
    save=false);

st = struct(inst=["0", "2"], tag=["Smoking status", "Smoking status MRI"]);
pth = fullfile(fileparts(which("phenoParser.m")), "UKB_PHENO"); % default path to save phenotypes

for k = 1:numel(st.tag)
    idx = tmp.tag.endsWith(st.inst(k));
    UKB_STRUCT_ALL = table2struct(tmp(idx, :));

    rmidx = UKB_STRUCT_ALL.rawUKB < 0; % prefer not to answer

    UKB_STRUCT_ALL.exeid = UKB_STRUCT_ALL.eid(rmidx);
    UKB_STRUCT_ALL.rawUKB(rmidx) = [];
    UKB_STRUCT_ALL.eid(rmidx) = [];
    UKB_STRUCT_ALL.termMeaning(rmidx) = [];
    UKB_STRUCT_ALL.tag = st.tag(k);
    
    out_name = fullfile(pth, matlab.lang.makeValidName(st.tag(k)) + ".mat");
    save(out_name, "UKB_STRUCT_ALL")
    clear UKB_STRUCT_ALL

end


end % END

%% subfunctions ===========================================================
function defDiabetes

codes = join(["E11" + (0:9), "E14" + (0:9)], ";");
df1 = phenoParser(query=["1220;1223", ...
    codes, "3", ...
    "1140868902;1140857584;1141171652;1141156984;1141189094;1141177606;1140874740;1140874706;1140874724;1140874736;1140874712;1140874746;1140857586;1140874728;1140874718;1140874650;1140857494;1140874744;1141152590;1140874646;1141157284;1140874658;1140921964;1140868908;1140874686;1140874660;1140857496;1140910564;1140910566;1140874678;1140874716;1140857500;1140883066;1140857590;1140874732;1140884600;1140874652;1141173882;1141168668;1140874690;1140882964;1141171646;1140874680;1141168660;1141177600;1141189090;1141173786;1140874666;1140874664;1140874674;1141153254;1141171508"],...
    tag=["df20002", "icd10", "medication", "drugs20003"], ...
    df=["20002", "", "6153;6177", "20003"], ...
    coding=["", "19", "", ""], ...
    save=false,...
    surv=false, ...
    all=true, ...
    merge=true);

df2 = phenoParser(query="",...
    df="30740;30750", ...
    instance="0", ...
    save=false,...
    surv=false, ...
    all=true);


% HbA1c ≥ 48 mmol/mol (6.5%) and Glucose ≥ 11.1 mmol/L (200 mg/dL)
df2.source = cell(height(df2), 1);
for k = 1:height(df2)
    if df2.tag(k).lower.contains("glucose")
        idx = df2.rawUKB{k} >= 11.1;
        rawtxt = "glucose>=11.1 mmol/L";
    elseif df2.tag(k).lower.contains("hba1c")
        idx = df2.rawUKB{k} >= 48;
        rawtxt = "HbA1c ≥ 48 mmol/mol";
    else
        error("only Glucose/HbA1c should be read!")
    end

    df2.eid{k}(~idx) = [];
    [df2.rawUKB{k}, df2.termMeaning{k}] = deal(repmat(rawtxt, nnz(idx), 1));
    df2.source{k} = repmat(df2.info(k).dfraw, nnz(idx), 1);
    df2.numericFlag(k) = false;
end
df2.unit = [];

if isstruct(df1)
    df1 = struct2table(df1, AsArray=true);
end

% create a composite df for all data-fields used to define type 2 diabetes
df_all1 = df1.info.df;
df_all2 = arrayfun(@(x)x.dfraw, df2.info);
df_all = df_all1 + "," + df_all2.join(",");
df_date = df2.info(1).date;

df1.info = []; df2.info = [];
df = [df1; df2]; clear df1 df2

% merge all chunks into a single pheno struct
UKB_STRUCT_ALL = struct;
UKB_STRUCT_ALL.tag = "Diabetes";
UKB_STRUCT_ALL.numericFlag = false;
UKB_STRUCT_ALL.eid = unique(vertcat(df.eid{:}));
[UKB_STRUCT_ALL.rawUKB, UKB_STRUCT_ALL.termMeaning, ...
    UKB_STRUCT_ALL.source] = deal(strings(numel(UKB_STRUCT_ALL.eid), 1));
UKB_STRUCT_ALL.info.date = df_date;
UKB_STRUCT_ALL.info.basket = "RAP";
UKB_STRUCT_ALL.info.df = df_all;

for j = 1:height(df)    
    [idx1, idx2] = ismember(UKB_STRUCT_ALL.eid, df.eid{j}); idx2(idx2 < 1) = [];
    fi = ["rawUKB", "source", "termMeaning"];
    for k = 1:numel(fi)
        UKB_STRUCT_ALL.(fi(k))(idx1) = df.(fi(k)){j}(idx2) + "," + UKB_STRUCT_ALL.(fi(k))(idx1);
    end
    
end

UKB_STRUCT_ALL.rawUKB = regexprep(UKB_STRUCT_ALL.rawUKB, '(^,$|,$|^,)', '');
UKB_STRUCT_ALL.termMeaning = regexprep(UKB_STRUCT_ALL.termMeaning, '(^,$|,$|^,)', '');
UKB_STRUCT_ALL.source = regexprep(UKB_STRUCT_ALL.source, '(^,$|,$|^,)', '');

% default path to save phenotypes
pth = fullfile(fileparts(which("phenoParser.m")), "UKB_PHENO");
save(fullfile(pth, "Diabetes.mat"), "UKB_STRUCT_ALL")

end % END
%%%%%%%%%%%%%%%%%   ALOFT   %%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all

clearvars -except D1* D2* sig_etv*

start=pwd;

num_mems=40; %for indexing
num_times=5;

spatial_agg = 0 %1 for yes, 0 for no
outofpoly=1 %1 does outside of poly, 0 inside

% poly_y = [28,42,42,28]; %sw,nw,ne,se corners 
% poly_x = [-106,-106,-100,-100]; %sw,nw,ne,se corners

%domain 2 LBCS (200 km)
poly_y = [30.24, 39.67, 39.7, 30.27]; %sw,nw,ne,se corners 
poly_x = [-101.73, -102.37, -93.9, -94.56]; %sw,nw,ne,se corners


% D1varD2var_etvheight_T = etvheight_T;
% D1varD2var_etvheight_TD = etvheight_TD;
% D1varD2var_etvheight_WIND = etvheight_WIND;
% D1varD2var_rmseheight_T = rmseheight_T;
% D1varD2var_rmseheight_TD = rmseheight_TD;
% D1varD2var_rmseheight_WIND = rmseheight_WIND;
% D1varD2var_varianceheight_T = varianceheight_T;
% D1varD2var_varianceheight_TD = varianceheight_TD;
% D1varD2var_varianceheight_WIND = varianceheight_WIND;
% D1varD2var_biasheight_T = biasheight_T;
% D1varD2var_biasheight_TD = biasheight_TD;
% D1varD2var_biasheight_WIND = biasheight_WIND;


% D2PBL_etvheight_T = etvheight_T;
% D2PBL_etvheight_TD = etvheight_TD;
% D2PBL_etvheight_WIND = etvheight_WIND;
% D2PBL_rmseheight_T = rmseheight_T;
% D2PBL_rmseheight_TD = rmseheight_TD;
% D2PBL_rmseheight_WIND = rmseheight_WIND;
% D2PBL_varianceheight_T = varianceheight_T;
% D2PBL_varianceheight_TD = varianceheight_TD;
% D2PBL_varianceheight_WIND = varianceheight_WIND;
% D2PBL_biasheight_T = biasheight_T;
% D2PBL_biasheight_TD = biasheight_TD;
% D2PBL_biasheight_WIND = biasheight_WIND;


% we=2.5; %wind error sigma 2.5 m/s?
%%%%%%%%%%%%%%%%%%% OBSERVATION ERROR ST DEV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_TD150=0.9; sigma_T150=0.9; sigma_U150=3.0; sigma_V150=sigma_U150; sigma_WIND150=sigma_U150;
sigma_TD200=0.9; sigma_T200=0.9; sigma_U200=3.0; sigma_V200=sigma_U200; sigma_WIND200=sigma_U200;
sigma_TD250=0.9; sigma_T250=0.9; sigma_U250=3.0; sigma_V250=sigma_U250; sigma_WIND250=sigma_U250;
sigma_TD300=0.8; sigma_T300=0.9; sigma_U300=3.0; sigma_V300=sigma_U300; sigma_WIND300=sigma_U300;
sigma_TD400=0.8; sigma_T400=0.8; sigma_U400=2.6; sigma_V400=sigma_U400; sigma_WIND400=sigma_U400;
sigma_TD500=0.8; sigma_T500=0.8; sigma_U500=2.1; sigma_V500=sigma_U500; sigma_WIND500=sigma_U500;
sigma_TD700=0.8; sigma_T700=0.8; sigma_U700=1.6; sigma_V700=sigma_U700; sigma_WIND700=sigma_U700;
sigma_TD850=0.8; sigma_T850=0.8; sigma_U850=1.5; sigma_V850=sigma_U850; sigma_WIND850=sigma_U850;
sigma_TD925=1.0; sigma_T925=1.0; sigma_U925=1.5; sigma_V925=sigma_U925; sigma_WIND925=sigma_U925;
sigma_TD1000=1.0; sigma_T1000=1.0; sigma_U1000=1.5; sigma_V1000=sigma_U1000; sigma_WIND1000=sigma_U1000;

%%%%%%%%%%%%%%%%%%%% ALLOCATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rank_sum_TD150=zeros(num_mems+1,num_times); rank_sum_T150=zeros(num_mems+1,num_times); rank_sum_WIND150=zeros(num_mems+1,num_times);
rank_sum_TD200=zeros(num_mems+1,num_times); rank_sum_T200=zeros(num_mems+1,num_times); rank_sum_WIND200=zeros(num_mems+1,num_times);
rank_sum_TD250=zeros(num_mems+1,num_times); rank_sum_T250=zeros(num_mems+1,num_times); rank_sum_WIND250=zeros(num_mems+1,num_times);
rank_sum_TD300=zeros(num_mems+1,num_times); rank_sum_T300=zeros(num_mems+1,num_times); rank_sum_WIND300=zeros(num_mems+1,num_times);
rank_sum_TD400=zeros(num_mems+1,num_times); rank_sum_T400=zeros(num_mems+1,num_times); rank_sum_WIND400=zeros(num_mems+1,num_times);
rank_sum_TD500=zeros(num_mems+1,num_times); rank_sum_T500=zeros(num_mems+1,num_times); rank_sum_WIND500=zeros(num_mems+1,num_times);
rank_sum_TD700=zeros(num_mems+1,num_times); rank_sum_T700=zeros(num_mems+1,num_times); rank_sum_WIND700=zeros(num_mems+1,num_times);
rank_sum_TD850=zeros(num_mems+1,num_times); rank_sum_T850=zeros(num_mems+1,num_times); rank_sum_WIND850=zeros(num_mems+1,num_times);
rank_sum_TD925=zeros(num_mems+1,num_times); rank_sum_T925=zeros(num_mems+1,num_times); rank_sum_WIND925=zeros(num_mems+1,num_times);
rank_sum_TD1000=zeros(num_mems+1,num_times); rank_sum_T1000=zeros(num_mems+1,num_times); rank_sum_WIND1000=zeros(num_mems+1,num_times);

type='var'; %either 'fixed','variedD1D2' or 'diff'
typed1='fix';
typed1_dir='PBLvar';

domain_list=[
%         'D1'
        'D2'
        ]; %only do 1 domain at a time
list=[
    '2015050600'
    '2015050612'
    '2015050700'
    '2015050712'
    '2015050800'
    '2015050812'
    '2015050900'
    '2015050912'   
    ];

T150=[]; TD150=[]; U150=[]; V150=[];
T200=[]; TD200=[]; U200=[]; V200=[];
T250=[]; TD250=[]; U250=[]; V250=[];
T300=[]; TD300=[]; U300=[]; V300=[];
T400=[]; TD400=[]; U400=[]; V400=[];
T500=[]; TD500=[]; U500=[]; V500=[];
T700=[]; TD700=[]; U700=[]; V700=[];
T850=[]; TD850=[]; U850=[]; V850=[];
T925=[]; TD925=[]; U925=[]; V925=[];
T1000=[]; TD1000=[]; U1000=[]; V1000=[];

latT150=[]; lonT150=[]; latTD150=[]; lonTD150=[]; latU150=[]; lonU150=[];
latT200=[]; lonT200=[]; latTD200=[]; lonTD200=[]; latU200=[]; lonU200=[];
latT250=[]; lonT250=[]; latTD250=[]; lonTD250=[]; latU250=[]; lonU250=[];
latT300=[]; lonT300=[]; latTD300=[]; lonTD300=[]; latU300=[]; lonU300=[];
latT400=[]; lonT400=[]; latTD400=[]; lonTD400=[]; latU400=[]; lonU400=[];
latT500=[]; lonT500=[]; latTD500=[]; lonTD500=[]; latU500=[]; lonU500=[];
latT700=[]; lonT700=[]; latTD700=[]; lonTD700=[]; latU700=[]; lonU700=[];
latT850=[]; lonT850=[]; latTD850=[]; lonTD850=[]; latU850=[]; lonU850=[];
latT925=[]; lonT925=[]; latTD925=[]; lonTD925=[]; latU925=[]; lonU925=[];
latT1000=[]; lonT1000=[]; latTD1000=[]; lonTD1000=[]; latU1000=[]; lonU1000=[];


% cd('/Users/brock/Google Drive/verif_files/D1var')

for y=1:size(domain_list,1)
    domain=domain_list(y,:);
cd([start '/D1' typed1_dir])
display(['Working in directory: ' typed1_dir])    
    
for z=1:size(list,1);

datestring=list(z,:);
% cd(['./' list(z,1:8) '_' dir_type 'physics'])

for t=1:num_times;
    
    fname=['verifaloft' domain type '_' datestring '_f' num2str((t-1)*12)] 
    ncid=netcdf.open(fname,'NOWRITE');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMPERATURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%% Have to preserve full members-obs for rank histograms %%%


    varid=netcdf.inqVarID(ncid,'LAT');
    LAT=netcdf.getVar(ncid,varid,'double');
    
    varid=netcdf.inqVarID(ncid,'LON');
    LON=netcdf.getVar(ncid,varid,'double');



    %%% T150
    varid=netcdf.inqVarID(ncid,'T150');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatT150=LAT;    blahlatT150(blahinds,:)=[];
    blahlonT150=LON;    blahlonT150(blahinds,:)=[];
    latT150(size(latT150,1)+1:size(latT150,1)+size(blah,1),1) = blahlatT150;
    lonT150(size(lonT150,1)+1:size(lonT150,1)+size(blah,1),1) = blahlonT150;    
    T150(size(T150,1)+1:size(T150,1)+size(blah,1),:,1) = blah;

    %%% T200
    varid=netcdf.inqVarID(ncid,'T200');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatT200=LAT;    blahlatT200(blahinds,:)=[];
    blahlonT200=LON;    blahlonT200(blahinds,:)=[];
    latT200(size(latT200,1)+1:size(latT200,1)+size(blah,1),1) = blahlatT200;
    lonT200(size(lonT200,1)+1:size(lonT200,1)+size(blah,1),1) = blahlonT200;    
    T200(size(T200,1)+1:size(T200,1)+size(blah,1),:,1) = blah;
    
    %%% T250
    varid=netcdf.inqVarID(ncid,'T250');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatT250=LAT;    blahlatT250(blahinds,:)=[];
    blahlonT250=LON;    blahlonT250(blahinds,:)=[];
    latT250(size(latT250,1)+1:size(latT250,1)+size(blah,1),1) = blahlatT250;
    lonT250(size(lonT250,1)+1:size(lonT250,1)+size(blah,1),1) = blahlonT250;    
    T250(size(T250,1)+1:size(T250,1)+size(blah,1),:,1) = blah;
    
    %%% T300
    varid=netcdf.inqVarID(ncid,'T300');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatT300=LAT;    blahlatT300(blahinds,:)=[];
    blahlonT300=LON;    blahlonT300(blahinds,:)=[];
    latT300(size(latT300,1)+1:size(latT300,1)+size(blah,1),1) = blahlatT300;
    lonT300(size(lonT300,1)+1:size(lonT300,1)+size(blah,1),1) = blahlonT300;    
    T300(size(T300,1)+1:size(T300,1)+size(blah,1),:,1) = blah;  
    
    %%% T400
    varid=netcdf.inqVarID(ncid,'T400');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatT400=LAT;    blahlatT400(blahinds,:)=[];
    blahlonT400=LON;    blahlonT400(blahinds,:)=[];
    latT400(size(latT400,1)+1:size(latT400,1)+size(blah,1),1) = blahlatT400;
    lonT400(size(lonT400,1)+1:size(lonT400,1)+size(blah,1),1) = blahlonT400;    
    T400(size(T400,1)+1:size(T400,1)+size(blah,1),:,1) = blah;  
    
    %%% T500
    varid=netcdf.inqVarID(ncid,'T500');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatT500=LAT;    blahlatT500(blahinds,:)=[];
    blahlonT500=LON;    blahlonT500(blahinds,:)=[];
    latT500(size(latT500,1)+1:size(latT500,1)+size(blah,1),1) = blahlatT500;
    lonT500(size(lonT500,1)+1:size(lonT500,1)+size(blah,1),1) = blahlonT500;    
    T500(size(T500,1)+1:size(T500,1)+size(blah,1),:,1) = blah;       
    
    
    %%% T700
    varid=netcdf.inqVarID(ncid,'T700');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatT700=LAT;    blahlatT700(blahinds,:)=[];
    blahlonT700=LON;    blahlonT700(blahinds,:)=[];
    latT700(size(latT700,1)+1:size(latT700,1)+size(blah,1),1) = blahlatT700;
    lonT700(size(lonT700,1)+1:size(lonT700,1)+size(blah,1),1) = blahlonT700;    
    T700(size(T700,1)+1:size(T700,1)+size(blah,1),:,1) = blah;  
    
    %%% T850
    varid=netcdf.inqVarID(ncid,'T850');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatT850=LAT;    blahlatT850(blahinds,:)=[];
    blahlonT850=LON;    blahlonT850(blahinds,:)=[];
    latT850(size(latT850,1)+1:size(latT850,1)+size(blah,1),1) = blahlatT850;
    lonT850(size(lonT850,1)+1:size(lonT850,1)+size(blah,1),1) = blahlonT850;    
    T850(size(T850,1)+1:size(T850,1)+size(blah,1),:,1) = blah;       
    
    %%% T925
    varid=netcdf.inqVarID(ncid,'T925');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatT925=LAT;    blahlatT925(blahinds,:)=[];
    blahlonT925=LON;    blahlonT925(blahinds,:)=[];
    latT925(size(latT925,1)+1:size(latT925,1)+size(blah,1),1) = blahlatT925;
    lonT925(size(lonT925,1)+1:size(lonT925,1)+size(blah,1),1) = blahlonT925;    
    T925(size(T925,1)+1:size(T925,1)+size(blah,1),:,1) = blah;   
    
    %%% T1000
    varid=netcdf.inqVarID(ncid,'T1000');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatT1000=LAT;    blahlatT1000(blahinds,:)=[];
    blahlonT1000=LON;    blahlonT1000(blahinds,:)=[];
    latT1000(size(latT1000,1)+1:size(latT1000,1)+size(blah,1),1) = blahlatT1000;
    lonT1000(size(lonT1000,1)+1:size(lonT1000,1)+size(blah,1),1) = blahlonT1000;    
    T1000(size(T1000,1)+1:size(T1000,1)+size(blah,1),:,1) = blah;       
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%%%%%%%%%%%%%%%%% DEWPOINT TEMPERATURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% TD150
    varid=netcdf.inqVarID(ncid,'TD150');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatTD150=LAT;    blahlatTD150(blahinds,:)=[];
    blahlonTD150=LON;    blahlonTD150(blahinds,:)=[];
    latTD150(size(latTD150,1)+1:size(latTD150,1)+size(blah,1),1) = blahlatTD150;
    lonTD150(size(lonTD150,1)+1:size(lonTD150,1)+size(blah,1),1) = blahlonTD150;    
    TD150(size(TD150,1)+1:size(TD150,1)+size(blah,1),:,1) = blah;

    %%% TD200
    varid=netcdf.inqVarID(ncid,'TD200');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatTD200=LAT;    blahlatTD200(blahinds,:)=[];
    blahlonTD200=LON;    blahlonTD200(blahinds,:)=[];
    latTD200(size(latTD200,1)+1:size(latTD200,1)+size(blah,1),1) = blahlatTD200;
    lonTD200(size(lonTD200,1)+1:size(lonTD200,1)+size(blah,1),1) = blahlonTD200;    
    TD200(size(TD200,1)+1:size(TD200,1)+size(blah,1),:,1) = blah;
    
    %%% TD250
    varid=netcdf.inqVarID(ncid,'TD250');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatTD250=LAT;    blahlatTD250(blahinds,:)=[];
    blahlonTD250=LON;    blahlonTD250(blahinds,:)=[];
    latTD250(size(latTD250,1)+1:size(latTD250,1)+size(blah,1),1) = blahlatTD250;
    lonTD250(size(lonTD250,1)+1:size(lonTD250,1)+size(blah,1),1) = blahlonTD250;    
    TD250(size(TD250,1)+1:size(TD250,1)+size(blah,1),:,1) = blah;
    
    %%% TD300
    varid=netcdf.inqVarID(ncid,'TD300');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatTD300=LAT;    blahlatTD300(blahinds,:)=[];
    blahlonTD300=LON;    blahlonTD300(blahinds,:)=[];
    latTD300(size(latTD300,1)+1:size(latTD300,1)+size(blah,1),1) = blahlatTD300;
    lonTD300(size(lonTD300,1)+1:size(lonTD300,1)+size(blah,1),1) = blahlonTD300;    
    TD300(size(TD300,1)+1:size(TD300,1)+size(blah,1),:,1) = blah;  
    
    %%% TD400
    varid=netcdf.inqVarID(ncid,'TD400');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatTD400=LAT;    blahlatTD400(blahinds,:)=[];
    blahlonTD400=LON;    blahlonTD400(blahinds,:)=[];
    latTD400(size(latTD400,1)+1:size(latTD400,1)+size(blah,1),1) = blahlatTD400;
    lonTD400(size(lonTD400,1)+1:size(lonTD400,1)+size(blah,1),1) = blahlonTD400;    
    TD400(size(TD400,1)+1:size(TD400,1)+size(blah,1),:,1) = blah;  
    
    %%% TD500
    varid=netcdf.inqVarID(ncid,'TD500');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatTD500=LAT;    blahlatTD500(blahinds,:)=[];
    blahlonTD500=LON;    blahlonTD500(blahinds,:)=[];
    latTD500(size(latTD500,1)+1:size(latTD500,1)+size(blah,1),1) = blahlatTD500;
    lonTD500(size(lonTD500,1)+1:size(lonTD500,1)+size(blah,1),1) = blahlonTD500;    
    TD500(size(TD500,1)+1:size(TD500,1)+size(blah,1),:,1) = blah;       
    
    %%% TD700
    varid=netcdf.inqVarID(ncid,'TD700');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatTD700=LAT;    blahlatTD700(blahinds,:)=[];
    blahlonTD700=LON;    blahlonTD700(blahinds,:)=[];
    latTD700(size(latTD700,1)+1:size(latTD700,1)+size(blah,1),1) = blahlatTD700;
    lonTD700(size(lonTD700,1)+1:size(lonTD700,1)+size(blah,1),1) = blahlonTD700;    
    TD700(size(TD700,1)+1:size(TD700,1)+size(blah,1),:,1) = blah;  
    
    %%% TD850
    varid=netcdf.inqVarID(ncid,'TD850');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatTD850=LAT;    blahlatTD850(blahinds,:)=[];
    blahlonTD850=LON;    blahlonTD850(blahinds,:)=[];
    latTD850(size(latTD850,1)+1:size(latTD850,1)+size(blah,1),1) = blahlatTD850;
    lonTD850(size(lonTD850,1)+1:size(lonTD850,1)+size(blah,1),1) = blahlonTD850;    
    TD850(size(TD850,1)+1:size(TD850,1)+size(blah,1),:,1) = blah;       
    
    %%% TD925
    varid=netcdf.inqVarID(ncid,'TD925');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatTD925=LAT;    blahlatTD925(blahinds,:)=[];
    blahlonTD925=LON;    blahlonTD925(blahinds,:)=[];
    latTD925(size(latTD925,1)+1:size(latTD925,1)+size(blah,1),1) = blahlatTD925;
    lonTD925(size(lonTD925,1)+1:size(lonTD925,1)+size(blah,1),1) = blahlonTD925;    
    TD925(size(TD925,1)+1:size(TD925,1)+size(blah,1),:,1) = blah;   
    
    %%% TD1000
    varid=netcdf.inqVarID(ncid,'TD1000');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<=0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatTD1000=LAT;    blahlatTD1000(blahinds,:)=[];
    blahlonTD1000=LON;    blahlonTD1000(blahinds,:)=[];
    latTD1000(size(latTD1000,1)+1:size(latTD1000,1)+size(blah,1),1) = blahlatTD1000;
    lonTD1000(size(lonTD1000,1)+1:size(lonTD1000,1)+size(blah,1),1) = blahlonTD1000;    
    TD1000(size(TD1000,1)+1:size(TD1000,1)+size(blah,1),:,1) = blah;       




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% U-WIND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    %%% U150
    varid=netcdf.inqVarID(ncid,'U150');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatU150=LAT;    blahlatU150(blahinds,:)=[];
    blahlonU150=LON;    blahlonU150(blahinds,:)=[];
    latU150(size(latU150,1)+1:size(latU150,1)+size(blah,1),1) = blahlatU150;
    lonU150(size(lonU150,1)+1:size(lonU150,1)+size(blah,1),1) = blahlonU150;    
    U150(size(U150,1)+1:size(U150,1)+size(blah,1),:,1) = blah;

    %%% U200
    varid=netcdf.inqVarID(ncid,'U200');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatU200=LAT;    blahlatU200(blahinds,:)=[];
    blahlonU200=LON;    blahlonU200(blahinds,:)=[];
    latU200(size(latU200,1)+1:size(latU200,1)+size(blah,1),1) = blahlatU200;
    lonU200(size(lonU200,1)+1:size(lonU200,1)+size(blah,1),1) = blahlonU200;    
    U200(size(U200,1)+1:size(U200,1)+size(blah,1),:,1) = blah;
    
    %%% U250
    varid=netcdf.inqVarID(ncid,'U250');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatU250=LAT;    blahlatU250(blahinds,:)=[];
    blahlonU250=LON;    blahlonU250(blahinds,:)=[];
    latU250(size(latU250,1)+1:size(latU250,1)+size(blah,1),1) = blahlatU250;
    lonU250(size(lonU250,1)+1:size(lonU250,1)+size(blah,1),1) = blahlonU250;    
    U250(size(U250,1)+1:size(U250,1)+size(blah,1),:,1) = blah;
    
    %%% U300
    varid=netcdf.inqVarID(ncid,'U300');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatU300=LAT;    blahlatU300(blahinds,:)=[];
    blahlonU300=LON;    blahlonU300(blahinds,:)=[];
    latU300(size(latU300,1)+1:size(latU300,1)+size(blah,1),1) = blahlatU300;
    lonU300(size(lonU300,1)+1:size(lonU300,1)+size(blah,1),1) = blahlonU300;    
    U300(size(U300,1)+1:size(U300,1)+size(blah,1),:,1) = blah;  
    
    %%% U400
    varid=netcdf.inqVarID(ncid,'U400');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatU400=LAT;    blahlatU400(blahinds,:)=[];
    blahlonU400=LON;    blahlonU400(blahinds,:)=[];
    latU400(size(latU400,1)+1:size(latU400,1)+size(blah,1),1) = blahlatU400;
    lonU400(size(lonU400,1)+1:size(lonU400,1)+size(blah,1),1) = blahlonU400;    
    U400(size(U400,1)+1:size(U400,1)+size(blah,1),:,1) = blah;  
    
    %%% U500
    varid=netcdf.inqVarID(ncid,'U500');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatU500=LAT;    blahlatU500(blahinds,:)=[];
    blahlonU500=LON;    blahlonU500(blahinds,:)=[];
    latU500(size(latU500,1)+1:size(latU500,1)+size(blah,1),1) = blahlatU500;
    lonU500(size(lonU500,1)+1:size(lonU500,1)+size(blah,1),1) = blahlonU500;    
    U500(size(U500,1)+1:size(U500,1)+size(blah,1),:,1) = blah;       
    
    
    %%% U700
    varid=netcdf.inqVarID(ncid,'U700');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatU700=LAT;    blahlatU700(blahinds,:)=[];
    blahlonU700=LON;    blahlonU700(blahinds,:)=[];
    latU700(size(latU700,1)+1:size(latU700,1)+size(blah,1),1) = blahlatU700;
    lonU700(size(lonU700,1)+1:size(lonU700,1)+size(blah,1),1) = blahlonU700;    
    U700(size(U700,1)+1:size(U700,1)+size(blah,1),:,1) = blah;  
    
    %%% U850
    varid=netcdf.inqVarID(ncid,'U850');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatU850=LAT;    blahlatU850(blahinds,:)=[];
    blahlonU850=LON;    blahlonU850(blahinds,:)=[];
    latU850(size(latU850,1)+1:size(latU850,1)+size(blah,1),1) = blahlatU850;
    lonU850(size(lonU850,1)+1:size(lonU850,1)+size(blah,1),1) = blahlonU850;    
    U850(size(U850,1)+1:size(U850,1)+size(blah,1),:,1) = blah;       
    
    %%% U925
    varid=netcdf.inqVarID(ncid,'U925');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatU925=LAT;    blahlatU925(blahinds,:)=[];
    blahlonU925=LON;    blahlonU925(blahinds,:)=[];
    latU925(size(latU925,1)+1:size(latU925,1)+size(blah,1),1) = blahlatU925;
    lonU925(size(lonU925,1)+1:size(lonU925,1)+size(blah,1),1) = blahlonU925;    
    U925(size(U925,1)+1:size(U925,1)+size(blah,1),:,1) = blah;   
    
    %%% U1000
    varid=netcdf.inqVarID(ncid,'U1000');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    blahlatU1000=LAT;    blahlatU1000(blahinds,:)=[];
    blahlonU1000=LON;    blahlonU1000(blahinds,:)=[];
    latU1000(size(latU1000,1)+1:size(latU1000,1)+size(blah,1),1) = blahlatU1000;
    lonU1000(size(lonU1000,1)+1:size(lonU1000,1)+size(blah,1),1) = blahlonU1000;    
    U1000(size(U1000,1)+1:size(U1000,1)+size(blah,1),:,1) = blah;      
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% V-WIND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    %%% V150
    varid=netcdf.inqVarID(ncid,'V150');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    latV150=LAT;    latV150(blahinds,:)=[];
    lonV150=LON;    lonV150(blahinds,:)=[];
    V150(size(V150,1)+1:size(V150,1)+size(blah,1),:,1) = blah;
%     V150(size(V150,1)+1:size(V150,1)+size(blah,1),:,2) = lat;
%     V150(size(V150,1)+1:size(V150,1)+size(blah,1),:,3) = lon;

    %%% V200
    varid=netcdf.inqVarID(ncid,'V200');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    V200(size(V200,1)+1:size(V200,1)+size(blah,1),:) = blah;
    latV200=LAT;    latV200(blahinds,:)=[];
    lonV200=LON;    lonV200(blahinds,:)=[];
    
    %%% V250
    varid=netcdf.inqVarID(ncid,'V250');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    V250(size(V250,1)+1:size(V250,1)+size(blah,1),:) = blah;    
    latV250=LAT;    latV250(blahinds,:)=[];
    lonV250=LON;    lonV250(blahinds,:)=[];
    
    %%% V300
    varid=netcdf.inqVarID(ncid,'V300');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    V300(size(V300,1)+1:size(V300,1)+size(blah,1),:) = blah;
    latV300=LAT;    latV300(blahinds,:)=[];
    lonV300=LON;    lonV300(blahinds,:)=[];  
    
    %%% V400
    varid=netcdf.inqVarID(ncid,'V400');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    V400(size(V400,1)+1:size(V400,1)+size(blah,1),:) = blah;
    latV400=LAT;    latV400(blahinds,:)=[];
    lonV400=LON;    lonV400(blahinds,:)=[];  
    
    %%% V500
    varid=netcdf.inqVarID(ncid,'V500');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    V500(size(V500,1)+1:size(V500,1)+size(blah,1),:) = blah;
    latV500=LAT;    latV500(blahinds,:)=[];
    lonV500=LON;    lonV500(blahinds,:)=[];       
    
    
    %%% V700
    varid=netcdf.inqVarID(ncid,'V700');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    V700(size(V700,1)+1:size(V700,1)+size(blah,1),:) = blah;
    latV700=LAT;    latV700(blahinds,:)=[];
    lonV700=LON;    lonV700(blahinds,:)=[];   
    
    %%% V850
    varid=netcdf.inqVarID(ncid,'V850');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    V850(size(V850,1)+1:size(V850,1)+size(blah,1),:) = blah;
    latV850=LAT;    latV850(blahinds,:)=[];
    lonV850=LON;    lonV850(blahinds,:)=[];       
    
    %%% V925
    varid=netcdf.inqVarID(ncid,'V925');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    V925(size(V925,1)+1:size(V925,1)+size(blah,1),:) = blah;
    latV925=LAT;    latV925(blahinds,:)=[];
    lonV925=LON;    lonV925(blahinds,:)=[];   
    
    %%% V1000
    varid=netcdf.inqVarID(ncid,'V1000');
    blah=netcdf.getVar(ncid,varid,'double');
    blah(blah<-9000 | blah==0)=nan;
    [row,col]=find(isnan(blah));
    blahinds=unique(row);
    blah(blahinds,:)=[];
    V1000(size(V1000,1)+1:size(V1000,1)+size(blah,1),:) = blah;
    latV1000=LAT;    latV1000(blahinds,:)=[];
    lonV1000=LON;    lonV1000(blahinds,:)=[]; 
      
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% DONE BUILDING ALL TIMES-CASES OBS VECTOR %%%
    clear ncid
end %times
end %cases
end %domain
    




%%% T150
modobsT150=zeros(size(T150,1),2); %Allocate
for i=1:size(T150,1)
modobsT150(i,1)=nanmean(T150(i,1:num_mems)); %take ens mean
modobsT150(i,2)=T150(i,num_mems+1); %no mean just raw obs value
end

%%% T200
modobsT200=zeros(size(T200,1),2); %Allocate
for i=1:size(T200,1)
modobsT200(i,1)=nanmean(T200(i,1:num_mems)); %take ens mean
modobsT200(i,2)=T200(i,num_mems+1); %no mean just raw obs value
end

%%% T250
modobsT250=zeros(size(T250,1),2); %Allocate
for i=1:size(T250,1)
modobsT250(i,1)=nanmean(T250(i,1:num_mems)); %take ens mean
modobsT250(i,2)=T250(i,num_mems+1); %no mean just raw obs value
end

%%% T300
modobsT300=zeros(size(T300,1),2); %Allocate
for i=1:size(T300,1)
modobsT300(i,1)=nanmean(T300(i,1:num_mems)); %take ens mean
modobsT300(i,2)=T300(i,num_mems+1); %no mean just raw obs value
end

%%% T400
modobsT400=zeros(size(T400,1),2); %Allocate
for i=1:size(T400,1)
modobsT400(i,1)=nanmean(T400(i,1:num_mems)); %take ens mean
modobsT400(i,2)=T400(i,num_mems+1); %no mean just raw obs value
end

%%% T500
modobsT500=zeros(size(T500,1),2); %Allocate
for i=1:size(T500,1)
modobsT500(i,1)=nanmean(T500(i,1:num_mems)); %take ens mean
modobsT500(i,2)=T500(i,num_mems+1); %no mean just raw obs value
end

%%% T700
modobsT700=zeros(size(T700,1),2); %Allocate
for i=1:size(T700,1)
modobsT700(i,1)=nanmean(T700(i,1:num_mems)); %take ens mean
modobsT700(i,2)=T700(i,num_mems+1); %no mean just raw obs value
end

%%% T850
modobsT850=zeros(size(T850,1),2); %Allocate
for i=1:size(T850,1)
modobsT850(i,1)=nanmean(T850(i,1:num_mems)); %take ens mean
modobsT850(i,2)=T850(i,num_mems+1); %no mean just raw obs value
end

%%% T925
modobsT925=zeros(size(T925,1),2); %Allocate
for i=1:size(T925,1)
modobsT925(i,1)=nanmean(T925(i,1:num_mems)); %take ens mean
modobsT925(i,2)=T925(i,num_mems+1); %no mean just raw obs value
end

%%% T1000
modobsT1000=zeros(size(T1000,1),2); %Allocate
for i=1:size(T1000,1)
modobsT1000(i,1)=nanmean(T1000(i,1:num_mems)); %take ens mean
modobsT1000(i,2)=T1000(i,num_mems+1); %no mean just raw obs value
end




%%% TD150
modobsTD150=zeros(size(TD150,1),2); %Allocate
for i=1:size(TD150,1)
modobsTD150(i,1)=nanmean(TD150(i,1:num_mems)); %take ens mean
modobsTD150(i,2)=TD150(i,num_mems+1); %no mean just raw obs value
end

%%% TD200
modobsTD200=zeros(size(TD200,1),2); %Allocate
for i=1:size(TD200,1)
modobsTD200(i,1)=nanmean(TD200(i,1:num_mems)); %take ens mean
modobsTD200(i,2)=TD200(i,num_mems+1); %no mean just raw obs value
end

%%% TD250
modobsTD250=zeros(size(TD250,1),2); %Allocate
for i=1:size(TD250,1)
modobsTD250(i,1)=nanmean(TD250(i,1:num_mems)); %take ens mean
modobsTD250(i,2)=TD250(i,num_mems+1); %no mean just raw obs value
end

%%% TD300
modobsTD300=zeros(size(TD300,1),2); %Allocate
for i=1:size(TD300,1)
modobsTD300(i,1)=nanmean(TD300(i,1:num_mems)); %take ens mean
modobsTD300(i,2)=TD300(i,num_mems+1); %no mean just raw obs value
end

%%% TD400
modobsTD400=zeros(size(TD400,1),2); %Allocate
for i=1:size(TD400,1)
modobsTD400(i,1)=nanmean(TD400(i,1:num_mems)); %take ens mean
modobsTD400(i,2)=TD400(i,num_mems+1); %no mean just raw obs value
end

%%% TD500
modobsTD500=zeros(size(TD500,1),2); %Allocate
for i=1:size(TD500,1)
modobsTD500(i,1)=nanmean(TD500(i,1:num_mems)); %take ens mean
modobsTD500(i,2)=TD500(i,num_mems+1); %no mean just raw obs value
end

%%% TD700
modobsTD700=zeros(size(TD700,1),2); %Allocate
for i=1:size(TD700,1)
modobsTD700(i,1)=nanmean(TD700(i,1:num_mems)); %take ens mean
modobsTD700(i,2)=TD700(i,num_mems+1); %no mean just raw obs value
end

%%% TD850
modobsTD850=zeros(size(TD850,1),2); %Allocate
for i=1:size(TD850,1)
modobsTD850(i,1)=nanmean(TD850(i,1:num_mems)); %take ens mean
modobsTD850(i,2)=TD850(i,num_mems+1); %no mean just raw obs value
end

%%% TD925
modobsTD925=zeros(size(TD925,1),2); %Allocate
for i=1:size(TD925,1)
modobsTD925(i,1)=nanmean(TD925(i,1:num_mems)); %take ens mean
modobsTD925(i,2)=TD925(i,num_mems+1); %no mean just raw obs value
end

%%% TD1000
modobsTD1000=zeros(size(TD1000,1),2); %Allocate
for i=1:size(TD1000,1)
modobsTD1000(i,1)=nanmean(TD1000(i,1:num_mems)); %take ens mean
modobsTD1000(i,2)=TD1000(i,num_mems+1); %no mean just raw obs value
end



%%% U150
modobsU150=zeros(size(U150,1),2); %Allocate
for i=1:size(U150,1)
modobsU150(i,1)=nanmean(U150(i,1:num_mems)); %take ens mean
modobsU150(i,2)=U150(i,num_mems+1); %no mean just raw obs value
end

%%% U200
modobsU200=zeros(size(U200,1),2); %Allocate
for i=1:size(U200,1)
modobsU200(i,1)=nanmean(U200(i,1:num_mems)); %take ens mean
modobsU200(i,2)=U200(i,num_mems+1); %no mean just raw obs value
end

%%% U250
modobsU250=zeros(size(U250,1),2); %Allocate
for i=1:size(U250,1)
modobsU250(i,1)=nanmean(U250(i,1:num_mems)); %take ens mean
modobsU250(i,2)=U250(i,num_mems+1); %no mean just raw obs value
end

%%% U300
modobsU300=zeros(size(U300,1),2); %Allocate
for i=1:size(U300,1)
modobsU300(i,1)=nanmean(U300(i,1:num_mems)); %take ens mean
modobsU300(i,2)=U300(i,num_mems+1); %no mean just raw obs value
end

%%% U400
modobsU400=zeros(size(U400,1),2); %Allocate
for i=1:size(U400,1)
modobsU400(i,1)=nanmean(U400(i,1:num_mems)); %take ens mean
modobsU400(i,2)=U400(i,num_mems+1); %no mean just raw obs value
end

%%% U500
modobsU500=zeros(size(U500,1),2); %Allocate
for i=1:size(U500,1)
modobsU500(i,1)=nanmean(U500(i,1:num_mems)); %take ens mean
modobsU500(i,2)=U500(i,num_mems+1); %no mean just raw obs value
end

%%% U700
modobsU700=zeros(size(U700,1),2); %Allocate
for i=1:size(U700,1)
modobsU700(i,1)=nanmean(U700(i,1:num_mems)); %take ens mean
modobsU700(i,2)=U700(i,num_mems+1); %no mean just raw obs value
end

%%% U850
modobsU850=zeros(size(U850,1),2); %Allocate
for i=1:size(U850,1)
modobsU850(i,1)=nanmean(U850(i,1:num_mems)); %take ens mean
modobsU850(i,2)=U850(i,num_mems+1); %no mean just raw obs value
end

%%% U925
modobsU925=zeros(size(U925,1),2); %Allocate
for i=1:size(U925,1)
modobsU925(i,1)=nanmean(U925(i,1:num_mems)); %take ens mean
modobsU925(i,2)=U925(i,num_mems+1); %no mean just raw obs value
end

%%% U1000
modobsU1000=zeros(size(U1000,1),2); %Allocate
for i=1:size(U1000,1)
modobsU1000(i,1)=nanmean(U1000(i,1:num_mems)); %take ens mean
modobsU1000(i,2)=U1000(i,num_mems+1); %no mean just raw obs value
end


%%% V150
modobsV150=zeros(size(V150,1),2); %Allocate
for i=1:size(V150,1)
modobsV150(i,1)=nanmean(V150(i,1:num_mems)); %take ens mean
modobsV150(i,2)=V150(i,num_mems+1); %no mean just raw obs value
end

%%% V200
modobsV200=zeros(size(V200,1),2); %Allocate
for i=1:size(V200,1)
modobsV200(i,1)=nanmean(V200(i,1:num_mems)); %take ens mean
modobsV200(i,2)=V200(i,num_mems+1); %no mean just raw obs value
end

%%% V250
modobsV250=zeros(size(V250,1),2); %Allocate
for i=1:size(V250,1)
modobsV250(i,1)=nanmean(V250(i,1:num_mems)); %take ens mean
modobsV250(i,2)=V250(i,num_mems+1); %no mean just raw obs value
end

%%% V300
modobsV300=zeros(size(V300,1),2); %Allocate
for i=1:size(V300,1)
modobsV300(i,1)=nanmean(V300(i,1:num_mems)); %take ens mean
modobsV300(i,2)=V300(i,num_mems+1); %no mean just raw obs value
end

%%% V400
modobsV400=zeros(size(V400,1),2); %Allocate
for i=1:size(V400,1)
modobsV400(i,1)=nanmean(V400(i,1:num_mems)); %take ens mean
modobsV400(i,2)=V400(i,num_mems+1); %no mean just raw obs value
end

%%% V500
modobsV500=zeros(size(V500,1),2); %Allocate
for i=1:size(V500,1)
modobsV500(i,1)=nanmean(V500(i,1:num_mems)); %take ens mean
modobsV500(i,2)=V500(i,num_mems+1); %no mean just raw obs value
end

%%% V700
modobsV700=zeros(size(V700,1),2); %Allocate
for i=1:size(V700,1)
modobsV700(i,1)=nanmean(V700(i,1:num_mems)); %take ens mean
modobsV700(i,2)=V700(i,num_mems+1); %no mean just raw obs value
end

%%% V850
modobsV850=zeros(size(V850,1),2); %Allocate
for i=1:size(V850,1)
modobsV850(i,1)=nanmean(V850(i,1:num_mems)); %take ens mean
modobsV850(i,2)=V850(i,num_mems+1); %no mean just raw obs value
end

%%% V925
modobsV925=zeros(size(V925,1),2); %Allocate
for i=1:size(V925,1)
modobsV925(i,1)=nanmean(V925(i,1:num_mems)); %take ens mean
modobsV925(i,2)=V925(i,num_mems+1); %no mean just raw obs value
end

%%% V1000
modobsV1000=zeros(size(V1000,1),2); %Allocate
for i=1:size(V1000,1)
modobsV1000(i,1)=nanmean(V1000(i,1:num_mems)); %take ens mean
modobsV1000(i,2)=V1000(i,num_mems+1); %no mean just raw obs value
end



%%%%%%%%%%%%%%%%%%%%%%%%%%% Full Wind Speed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%% WIND150
WIND150=sqrt(U150.^2 + V150.^2);
WIND150(WIND150==0)=nan;
WIND150(isnan(WIND150(:,num_mems+1)),:)=[];

modobsWIND150=sqrt(modobsU150.^2 + modobsV150.^2);
modobsWIND150(modobsWIND150==0)=nan;
modobsWIND150(isnan(modobsWIND150(:,2)),:)=[];

%%% WIND200
WIND200=sqrt(U200.^2 + V200.^2);
WIND200(WIND200==0)=nan;
WIND200(isnan(WIND200(:,num_mems+1)),:)=[];

modobsWIND200=sqrt(modobsU200.^2 + modobsV200.^2);
modobsWIND200(modobsWIND200==0)=nan;
modobsWIND200(isnan(modobsWIND200(:,2)),:)=[];

%%% WIND250
WIND250=sqrt(U250.^2 + V250.^2);
WIND250(WIND250==0)=nan;
WIND250(isnan(WIND250(:,num_mems+1)),:)=[];

modobsWIND250=sqrt(modobsU250.^2 + modobsV250.^2);
modobsWIND250(modobsWIND250==0)=nan;
modobsWIND250(isnan(modobsWIND250(:,2)),:)=[];

%%% WIND300
WIND300=sqrt(U300.^2 + V300.^2);
WIND300(WIND300==0)=nan;
WIND300(isnan(WIND300(:,num_mems+1)),:)=[];

modobsWIND300=sqrt(modobsU300.^2 + modobsV300.^2);
modobsWIND300(modobsWIND300==0)=nan;
modobsWIND300(isnan(modobsWIND300(:,2)),:)=[];

%%% WIND400
WIND400=sqrt(U400.^2 + V400.^2);
WIND400(WIND400==0)=nan;
WIND400(isnan(WIND400(:,num_mems+1)),:)=[];

modobsWIND400=sqrt(modobsU400.^2 + modobsV400.^2);
modobsWIND400(modobsWIND400==0)=nan;
modobsWIND400(isnan(modobsWIND400(:,2)),:)=[];

%%% WIND500
WIND500=sqrt(U500.^2 + V500.^2);
WIND500(WIND500==0)=nan;
WIND500(isnan(WIND500(:,num_mems+1)),:)=[];

modobsWIND500=sqrt(modobsU500.^2 + modobsV500.^2);
modobsWIND500(modobsWIND500==0)=nan;
modobsWIND500(isnan(modobsWIND500(:,2)),:)=[];

%%% WIND700
WIND700=sqrt(U700.^2 + V700.^2);
WIND700(WIND700==0)=nan;
WIND700(isnan(WIND700(:,num_mems+1)),:)=[];

modobsWIND700=sqrt(modobsU700.^2 + modobsV700.^2);
modobsWIND700(modobsWIND700==0)=nan;
modobsWIND700(isnan(modobsWIND700(:,2)),:)=[];

%%% WIND850
WIND850=sqrt(U850.^2 + V850.^2);
WIND850(WIND850==0)=nan;
WIND850(isnan(WIND850(:,num_mems+1)),:)=[];

modobsWIND850=sqrt(modobsU850.^2 + modobsV850.^2);
modobsWIND850(modobsWIND850==0)=nan;
modobsWIND850(isnan(modobsWIND850(:,2)),:)=[];

%%% WIND925
WIND925=sqrt(U925.^2 + V925.^2);
WIND925(WIND925==0)=nan;
WIND925(isnan(WIND925(:,num_mems+1)),:)=[];

modobsWIND925=sqrt(modobsU925.^2 + modobsV925.^2);
modobsWIND925(modobsWIND925==0)=nan;
modobsWIND925(isnan(modobsWIND925(:,2)),:)=[];

%%% WIND1000
WIND1000=sqrt(U1000.^2 + V1000.^2);
WIND1000(WIND1000==0)=nan;
WIND1000(isnan(WIND1000(:,num_mems+1)),:)=[];

modobsWIND1000=sqrt(modobsU1000.^2 + modobsV1000.^2);
modobsWIND1000(modobsWIND1000==0)=nan;
modobsWIND1000(isnan(modobsWIND1000(:,2)),:)=[];

%set WIND lat lons from U-component
latWIND150 = latU150;   lonWIND150 = lonU150;   
latWIND200 = latU200;   lonWIND200 = lonU200;
latWIND250 = latU250;   lonWIND250 = lonU250;
latWIND300 = latU300;   lonWIND300 = lonU300;
latWIND400 = latU400;   lonWIND400 = lonU400;
latWIND500 = latU500;   lonWIND500 = lonU500;
latWIND700 = latU700;   lonWIND700 = lonU700;
latWIND850 = latU850;   lonWIND850 = lonU850;
latWIND925 = latU925;   lonWIND925 = lonU925;
latWIND1000 = latU1000; lonWIND1000 = lonU1000;

    
%%%%%%%%% SPATIAL DISAGREGATION %%%%%%%%%%%%%%
if spatial_agg == 1
    display('Doing spatial aggregation based upon input polygon points')

    if outofpoly==1
        %%% only obtain verification OUTSIDE box (near LBC)
                %%% TEMPERATURE
        display('Doing spatial aggregation outside of polygon')         
        [in,on] = inpolygon(lonT150,latT150,poly_x,poly_y);
        T150 = T150(~in,:);  modobsT150 = modobsT150(~in,:);  latT150 = latT150(~in);  lonT150 = lonT150(~in);
        clearvars in on

        [in,on] = inpolygon(lonT200,latT200,poly_x,poly_y);
        T200 = T200(~in,:);  modobsT200 = modobsT200(~in,:);  latT200 = latT200(~in);  lonT200 = lonT200(~in);
        clearvars in on

        [in,on] = inpolygon(lonT250,latT250,poly_x,poly_y);
        T250 = T250(~in,:);  modobsT250 = modobsT250(~in,:);  latT250 = latT250(~in);  lonT250 = lonT250(~in);
        clearvars in on

        [in,on] = inpolygon(lonT300,latT300,poly_x,poly_y);
        T300 = T300(~in,:);  modobsT300 = modobsT300(~in,:);  latT300 = latT300(~in);  lonT300 = lonT300(~in);
        clearvars in on

        [in,on] = inpolygon(lonT400,latT400,poly_x,poly_y);
        T400 = T400(~in,:);  modobsT400 = modobsT400(~in,:);  latT400 = latT400(~in);  lonT400 = lonT400(~in);
        clearvars in on

        [in,on] = inpolygon(lonT500,latT500,poly_x,poly_y);
        T500 = T500(~in,:);  modobsT500 = modobsT500(~in,:);  latT500 = latT500(~in);  lonT500 = lonT500(~in);
        clearvars in on

        [in,on] = inpolygon(lonT700,latT700,poly_x,poly_y);
        T700 = T700(~in,:);  modobsT700 = modobsT700(~in,:);  latT700 = latT700(~in);  lonT700 = lonT700(~in);
        clearvars in on

        [in,on] = inpolygon(lonT850,latT850,poly_x,poly_y);
        T850 = T850(~in,:);  modobsT850 = modobsT850(~in,:);  latT850 = latT850(~in);  lonT850 = lonT850(~in);
        clearvars in on

        [in,on] = inpolygon(lonT925,latT925,poly_x,poly_y);
        T925 = T925(~in,:);  modobsT925 = modobsT925(~in,:);  latT925 = latT925(~in);  lonT925 = lonT925(~in);
        clearvars in on

        [in,on] = inpolygon(lonT1000,latT1000,poly_x,poly_y);
        T1000 = T1000(~in,:);  modobsT1000 = modobsT1000(~in,:);  latT1000 = latT1000(~in);  lonT1000 = lonT1000(~in);
        clearvars in on

        %%% DEWPOINT TEMPERATURE
        [in,on] = inpolygon(lonTD150,latTD150,poly_x,poly_y);
        TD150 = TD150(~in,:);  modobsTD150 = modobsTD150(~in,:);  latTD150 = latTD150(~in);  lonTD150 = lonTD150(~in);
        clearvars in on

        [in,on] = inpolygon(lonTD200,latTD200,poly_x,poly_y);
        TD200 = TD200(~in,:);  modobsTD200 = modobsTD200(~in,:);  latTD200 = latTD200(~in);  lonTD200 = lonTD200(~in);
        clearvars in on

        [in,on] = inpolygon(lonTD250,latTD250,poly_x,poly_y);
        TD250 = TD250(~in,:);  modobsTD250 = modobsTD250(~in,:);  latTD250 = latTD250(~in);  lonTD250 = lonTD250(~in);
        clearvars in on

        [in,on] = inpolygon(lonTD300,latTD300,poly_x,poly_y);
        TD300 = TD300(~in,:);  modobsTD300 = modobsTD300(~in,:);  latTD300 = latTD300(~in);  lonTD300 = lonTD300(~in);
        clearvars in on

        [in,on] = inpolygon(lonTD400,latTD400,poly_x,poly_y);
        TD400 = TD400(~in,:);  modobsTD400 = modobsTD400(~in,:);  latTD400 = latTD400(~in);  lonTD400 = lonTD400(~in);
        clearvars in on

        [in,on] = inpolygon(lonTD500,latTD500,poly_x,poly_y);
        TD500 = TD500(~in,:);  modobsTD500 = modobsTD500(~in,:);  latTD500 = latTD500(~in);  lonTD500 = lonTD500(~in);
        clearvars in on

        [in,on] = inpolygon(lonTD700,latTD700,poly_x,poly_y);
        TD700 = TD700(~in,:);  modobsTD700 = modobsTD700(~in,:);  latTD700 = latTD700(~in);  lonTD700 = lonTD700(~in);
        clearvars in on

        [in,on] = inpolygon(lonTD850,latTD850,poly_x,poly_y);
        TD850 = TD850(~in,:);  modobsTD850 = modobsTD850(~in,:);  latTD850 = latTD850(~in);  lonTD850 = lonTD850(~in);
        clearvars in on

        [in,on] = inpolygon(lonTD925,latTD925,poly_x,poly_y);
        TD925 = TD925(~in,:);  modobsTD925 = modobsTD925(~in,:);  latTD925 = latTD925(~in);  lonTD925 = lonTD925(~in);
        clearvars in on

        [in,on] = inpolygon(lonTD1000,latTD1000,poly_x,poly_y);
        TD1000 = TD1000(~in,:);  modobsTD1000 = modobsTD1000(~in,:);  latTD1000 = latTD1000(~in);  lonTD1000 = lonTD1000(~in);
        clearvars in on

        %%% WIND SPEED
        [in,on] = inpolygon(lonWIND150,latWIND150,poly_x,poly_y);
        WIND150 = WIND150(~in,:);  modobsWIND150 = modobsWIND150(~in,:);  latWIND150 = latWIND150(~in);  lonWIND150 = lonWIND150(~in);
        clearvars in on

        [in,on] = inpolygon(lonWIND200,latWIND200,poly_x,poly_y);
        WIND200 = WIND200(~in,:);  modobsWIND200 = modobsWIND200(~in,:);  latWIND200 = latWIND200(~in);  lonWIND200 = lonWIND200(~in);
        clearvars in on

        [in,on] = inpolygon(lonWIND250,latWIND250,poly_x,poly_y);
        WIND250 = WIND250(~in,:);  modobsWIND250 = modobsWIND250(~in,:);  latWIND250 = latWIND250(~in);  lonWIND250 = lonWIND250(~in);
        clearvars in on

        [in,on] = inpolygon(lonWIND300,latWIND300,poly_x,poly_y);
        WIND300 = WIND300(~in,:);  modobsWIND300 = modobsWIND300(~in,:);  latWIND300 = latWIND300(~in);  lonWIND300 = lonWIND300(~in);
        clearvars in on

        [in,on] = inpolygon(lonWIND400,latWIND400,poly_x,poly_y);
        WIND400 = WIND400(~in,:);  modobsWIND400 = modobsWIND400(~in,:);  latWIND400 = latWIND400(~in);  lonWIND400 = lonWIND400(~in);
        clearvars in on

        [in,on] = inpolygon(lonWIND500,latWIND500,poly_x,poly_y);
        WIND500 = WIND500(~in,:);  modobsWIND500 = modobsWIND500(~in,:);  latWIND500 = latWIND500(~in);  lonWIND500 = lonWIND500(~in);
        clearvars in on

        [in,on] = inpolygon(lonWIND700,latWIND700,poly_x,poly_y);
        WIND700 = WIND700(~in,:);  modobsWIND700 = modobsWIND700(~in,:);  latWIND700 = latWIND700(~in);  lonWIND700 = lonWIND700(~in);
        clearvars in on

        [in,on] = inpolygon(lonWIND850,latWIND850,poly_x,poly_y);
        WIND850 = WIND850(~in,:);  modobsWIND850 = modobsWIND850(~in,:);  latWIND850 = latWIND850(~in);  lonWIND850 = lonWIND850(~in);
        clearvars in on

        [in,on] = inpolygon(lonWIND925,latWIND925,poly_x,poly_y);
        WIND925 = WIND925(~in,:);  modobsWIND925 = modobsWIND925(~in,:);  latWIND925 = latWIND925(~in);  lonWIND925 = lonWIND925(~in);
        clearvars in on

        [in,on] = inpolygon(lonWIND1000,latWIND1000,poly_x,poly_y);
        WIND1000 = WIND1000(~in,:);  modobsWIND1000 = modobsWIND1000(~in,:);  latWIND1000 = latWIND1000(~in);  lonWIND1000 = lonWIND1000(~in);
        clearvars in on

 
    else
        
        %%%only obtain verification WITHIN box
        %%% TEMPERATURE
        display('Doing spatial aggregation inside of polygon')
        [in,on] = inpolygon(lonT150,latT150,poly_x,poly_y);
        T150 = T150(in,:);  modobsT150 = modobsT150(in,:);  latT150 = latT150(in);  lonT150 = lonT150(in);
        clearvars in on

        [in,on] = inpolygon(lonT200,latT200,poly_x,poly_y);
        T200 = T200(in,:);  modobsT200 = modobsT200(in,:);  latT200 = latT200(in);  lonT200 = lonT200(in);
        clearvars in on

        [in,on] = inpolygon(lonT250,latT250,poly_x,poly_y);
        T250 = T250(in,:);  modobsT250 = modobsT250(in,:);  latT250 = latT250(in);  lonT250 = lonT250(in);
        clearvars in on

        [in,on] = inpolygon(lonT300,latT300,poly_x,poly_y);
        T300 = T300(in,:);  modobsT300 = modobsT300(in,:);  latT300 = latT300(in);  lonT300 = lonT300(in);
        clearvars in on

        [in,on] = inpolygon(lonT400,latT400,poly_x,poly_y);
        T400 = T400(in,:);  modobsT400 = modobsT400(in,:);  latT400 = latT400(in);  lonT400 = lonT400(in);
        clearvars in on

        [in,on] = inpolygon(lonT500,latT500,poly_x,poly_y);
        T500 = T500(in,:);  modobsT500 = modobsT500(in,:);  latT500 = latT500(in);  lonT500 = lonT500(in);
        clearvars in on

        [in,on] = inpolygon(lonT700,latT700,poly_x,poly_y);
        T700 = T700(in,:);  modobsT700 = modobsT700(in,:);  latT700 = latT700(in);  lonT700 = lonT700(in);
        clearvars in on

        [in,on] = inpolygon(lonT850,latT850,poly_x,poly_y);
        T850 = T850(in,:);  modobsT850 = modobsT850(in,:);  latT850 = latT850(in);  lonT850 = lonT850(in);
        clearvars in on

        [in,on] = inpolygon(lonT925,latT925,poly_x,poly_y);
        T925 = T925(in,:);  modobsT925 = modobsT925(in,:);  latT925 = latT925(in);  lonT925 = lonT925(in);
        clearvars in on

        [in,on] = inpolygon(lonT1000,latT1000,poly_x,poly_y);
        T1000 = T1000(in,:);  modobsT1000 = modobsT1000(in,:);  latT1000 = latT1000(in);  lonT1000 = lonT1000(in);
        clearvars in on

        %%% DEWPOINT TEMPERATURE
        [in,on] = inpolygon(lonTD150,latTD150,poly_x,poly_y);
        TD150 = TD150(in,:);  modobsTD150 = modobsTD150(in,:);  latTD150 = latTD150(in);  lonTD150 = lonTD150(in);
        clearvars in on

        [in,on] = inpolygon(lonTD200,latTD200,poly_x,poly_y);
        TD200 = TD200(in,:);  modobsTD200 = modobsTD200(in,:);  latTD200 = latTD200(in);  lonTD200 = lonTD200(in);
        clearvars in on

        [in,on] = inpolygon(lonTD250,latTD250,poly_x,poly_y);
        TD250 = TD250(in,:);  modobsTD250 = modobsTD250(in,:);  latTD250 = latTD250(in);  lonTD250 = lonTD250(in);
        clearvars in on

        [in,on] = inpolygon(lonTD300,latTD300,poly_x,poly_y);
        TD300 = TD300(in,:);  modobsTD300 = modobsTD300(in,:);  latTD300 = latTD300(in);  lonTD300 = lonTD300(in);
        clearvars in on

        [in,on] = inpolygon(lonTD400,latTD400,poly_x,poly_y);
        TD400 = TD400(in,:);  modobsTD400 = modobsTD400(in,:);  latTD400 = latTD400(in);  lonTD400 = lonTD400(in);
        clearvars in on

        [in,on] = inpolygon(lonTD500,latTD500,poly_x,poly_y);
        TD500 = TD500(in,:);  modobsTD500 = modobsTD500(in,:);  latTD500 = latTD500(in);  lonTD500 = lonTD500(in);
        clearvars in on

        [in,on] = inpolygon(lonTD700,latTD700,poly_x,poly_y);
        TD700 = TD700(in,:);  modobsTD700 = modobsTD700(in,:);  latTD700 = latTD700(in);  lonTD700 = lonTD700(in);
        clearvars in on

        [in,on] = inpolygon(lonTD850,latTD850,poly_x,poly_y);
        TD850 = TD850(in,:);  modobsTD850 = modobsTD850(in,:);  latTD850 = latTD850(in);  lonTD850 = lonTD850(in);
        clearvars in on

        [in,on] = inpolygon(lonTD925,latTD925,poly_x,poly_y);
        TD925 = TD925(in,:);  modobsTD925 = modobsTD925(in,:);  latTD925 = latTD925(in);  lonTD925 = lonTD925(in);
        clearvars in on

        [in,on] = inpolygon(lonTD1000,latTD1000,poly_x,poly_y);
        TD1000 = TD1000(in,:);  modobsTD1000 = modobsTD1000(in,:);  latTD1000 = latTD1000(in);  lonTD1000 = lonTD1000(in);
        clearvars in on

        %%% WIND SPEED
        [in,on] = inpolygon(lonWIND150,latWIND150,poly_x,poly_y);
        WIND150 = WIND150(in,:);  modobsWIND150 = modobsWIND150(in,:);  latWIND150 = latWIND150(in);  lonWIND150 = lonWIND150(in);
        clearvars in on

        [in,on] = inpolygon(lonWIND200,latWIND200,poly_x,poly_y);
        WIND200 = WIND200(in,:);  modobsWIND200 = modobsWIND200(in,:);  latWIND200 = latWIND200(in);  lonWIND200 = lonWIND200(in);
        clearvars in on

        [in,on] = inpolygon(lonWIND250,latWIND250,poly_x,poly_y);
        WIND250 = WIND250(in,:);  modobsWIND250 = modobsWIND250(in,:);  latWIND250 = latWIND250(in);  lonWIND250 = lonWIND250(in);
        clearvars in on

        [in,on] = inpolygon(lonWIND300,latWIND300,poly_x,poly_y);
        WIND300 = WIND300(in,:);  modobsWIND300 = modobsWIND300(in,:);  latWIND300 = latWIND300(in);  lonWIND300 = lonWIND300(in);
        clearvars in on

        [in,on] = inpolygon(lonWIND400,latWIND400,poly_x,poly_y);
        WIND400 = WIND400(in,:);  modobsWIND400 = modobsWIND400(in,:);  latWIND400 = latWIND400(in);  lonWIND400 = lonWIND400(in);
        clearvars in on

        [in,on] = inpolygon(lonWIND500,latWIND500,poly_x,poly_y);
        WIND500 = WIND500(in,:);  modobsWIND500 = modobsWIND500(in,:);  latWIND500 = latWIND500(in);  lonWIND500 = lonWIND500(in);
        clearvars in on

        [in,on] = inpolygon(lonWIND700,latWIND700,poly_x,poly_y);
        WIND700 = WIND700(in,:);  modobsWIND700 = modobsWIND700(in,:);  latWIND700 = latWIND700(in);  lonWIND700 = lonWIND700(in);
        clearvars in on

        [in,on] = inpolygon(lonWIND850,latWIND850,poly_x,poly_y);
        WIND850 = WIND850(in,:);  modobsWIND850 = modobsWIND850(in,:);  latWIND850 = latWIND850(in);  lonWIND850 = lonWIND850(in);
        clearvars in on

        [in,on] = inpolygon(lonWIND925,latWIND925,poly_x,poly_y);
        WIND925 = WIND925(in,:);  modobsWIND925 = modobsWIND925(in,:);  latWIND925 = latWIND925(in);  lonWIND925 = lonWIND925(in);
        clearvars in on

        [in,on] = inpolygon(lonWIND1000,latWIND1000,poly_x,poly_y);
        WIND1000 = WIND1000(in,:);  modobsWIND1000 = modobsWIND1000(in,:);  latWIND1000 = latWIND1000(in);  lonWIND1000 = lonWIND1000(in);
        clearvars in on
    end
   
    
else
    display('Doing verification on entire domain')
end






    
    
%     obs_count_T150(t,z) = size(T150,1); obs_count_TD150(t,z) = size(TD150,1); obs_count_WIND150(t,z) = size(WIND150,1);
%     obs_count_T200(t,z) = size(T200,1); obs_count_TD200(t,z) = size(TD200,1); obs_count_WIND200(t,z) = size(WIND200,1);
%     obs_count_T250(t,z) = size(T250,1); obs_count_TD250(t,z) = size(TD250,1); obs_count_WIND250(t,z) = size(WIND250,1);
%     obs_count_T300(t,z) = size(T300,1); obs_count_TD300(t,z) = size(TD300,1); obs_count_WIND300(t,z) = size(WIND300,1);
%     obs_count_T400(t,z) = size(T400,1); obs_count_TD400(t,z) = size(TD400,1); obs_count_WIND400(t,z) = size(WIND400,1);
%     obs_count_T500(t,z) = size(T500,1); obs_count_TD500(t,z) = size(TD500,1); obs_count_WIND500(t,z) = size(WIND500,1);
%     obs_count_T700(t,z) = size(T700,1); obs_count_TD700(t,z) = size(TD700,1); obs_count_WIND700(t,z) = size(WIND700,1);
%     obs_count_T850(t,z) = size(T850,1); obs_count_TD850(t,z) = size(TD850,1); obs_count_WIND850(t,z) = size(WIND850,1);
%     obs_count_T925(t,z) = size(T925,1); obs_count_TD925(t,z) = size(TD925,1); obs_count_WIND925(t,z) = size(WIND925,1);    
%     obs_count_T1000(t,z) = size(T1000,1); obs_count_TD1000(t,z) = size(TD1000,1); obs_count_WIND1000(t,z) = size(WIND1000,1);
    
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%     netcdf.close(ncid)
%     clear ncid
  %end
  %end
  %end
%%%% READ IN ALL TIMES INTO ONE ARRAY FOR PROPER MEAN %%%%
  
%%%%%%%%%%%%%%%% STATISTICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% RMSE individual mems %%%
    for e=1:num_mems
        
%     RMSGPH300_timemem(t,e)=sqrt(mean((GPH300(:,e)-GPH300(:,num_mems+1)).^2));
%     RMSGPH400_timemem(t,e)=sqrt(mean((GPH400(:,e)-GPH400(:,num_mems+1)).^2));
%     RMSGPH500_timemem(t,e)=sqrt(mean((GPH500(:,e)-GPH500(:,num_mems+1)).^2));
%     RMSGPH700_timemem(t,e)=sqrt(mean((GPH700(:,e)-GPH700(:,num_mems+1)).^2));
%     RMSGPH850_timemem(t,e)=sqrt(mean((GPH850(:,e)-GPH850(:,num_mems+1)).^2));
%     RMSGPH925_timemem(t,e)=sqrt(mean((GPH925(:,e)-GPH925(:,num_mems+1)).^2));
     
    RMST150_timemem(e)=sqrt(mean((T150(:,e)-T150(:,num_mems+1)).^2));
    RMST200_timemem(e)=sqrt(mean((T200(:,e)-T200(:,num_mems+1)).^2));
    RMST250_timemem(e)=sqrt(mean((T250(:,e)-T250(:,num_mems+1)).^2));
    RMST300_timemem(e)=sqrt(mean((T300(:,e)-T300(:,num_mems+1)).^2));
    RMST400_timemem(e)=sqrt(mean((T400(:,e)-T400(:,num_mems+1)).^2));
    RMST500_timemem(e)=sqrt(mean((T500(:,e)-T500(:,num_mems+1)).^2));
    RMST700_timemem(e)=sqrt(mean((T700(:,e)-T700(:,num_mems+1)).^2));
    RMST850_timemem(e)=sqrt(mean((T850(:,e)-T850(:,num_mems+1)).^2));
    RMST925_timemem(e)=sqrt(mean((T925(:,e)-T925(:,num_mems+1)).^2));
    RMST1000_timemem(e)=sqrt(mean((T1000(:,e)-T1000(:,num_mems+1)).^2));
    
    RMSWIND150_timemem(e)=sqrt(mean((WIND150(:,e)-WIND150(:,num_mems+1)).^2));
    RMSWIND200_timemem(e)=sqrt(mean((WIND200(:,e)-WIND200(:,num_mems+1)).^2));
    RMSWIND250_timemem(e)=sqrt(mean((WIND250(:,e)-WIND250(:,num_mems+1)).^2));
    RMSWIND300_timemem(e)=sqrt(mean((WIND300(:,e)-WIND300(:,num_mems+1)).^2));
    RMSWIND400_timemem(e)=sqrt(mean((WIND400(:,e)-WIND400(:,num_mems+1)).^2));
    RMSWIND500_timemem(e)=sqrt(mean((WIND500(:,e)-WIND500(:,num_mems+1)).^2));
    RMSWIND700_timemem(e)=sqrt(mean((WIND700(:,e)-WIND700(:,num_mems+1)).^2));
    RMSWIND850_timemem(e)=sqrt(mean((WIND850(:,e)-WIND850(:,num_mems+1)).^2));
    RMSWIND925_timemem(e)=sqrt(mean((WIND925(:,e)-WIND925(:,num_mems+1)).^2));
    RMSWIND1000_timemem(e)=sqrt(mean((WIND1000(:,e)-WIND1000(:,num_mems+1)).^2));
    
    end


    %%% RMSE mean, Variance %%%
    biasT150 = mean(modobsT150(:,1)) - mean(modobsT150(:,2));
    biasT200 = mean(modobsT200(:,1)) - mean(modobsT200(:,2));
    biasT250 = mean(modobsT250(:,1)) - mean(modobsT250(:,2));
    biasT300 = mean(modobsT300(:,1)) - mean(modobsT300(:,2));
    biasT400 = mean(modobsT400(:,1)) - mean(modobsT400(:,2));
    biasT500 = mean(modobsT500(:,1)) - mean(modobsT500(:,2));
    biasT700 = mean(modobsT700(:,1)) - mean(modobsT700(:,2));
    biasT850 = mean(modobsT850(:,1)) - mean(modobsT850(:,2));
    biasT925 = mean(modobsT925(:,1)) - mean(modobsT925(:,2));
    biasT1000 = mean(modobsT1000(:,1)) - mean(modobsT1000(:,2));
    
    biasTD150 = mean(modobsTD150(:,1)) - mean(modobsTD150(:,2));
    biasTD200 = mean(modobsTD200(:,1)) - mean(modobsTD200(:,2));
    biasTD250 = mean(modobsTD250(:,1)) - mean(modobsTD250(:,2));
    biasTD300 = mean(modobsTD300(:,1)) - mean(modobsTD300(:,2));
    biasTD400 = mean(modobsTD400(:,1)) - mean(modobsTD400(:,2));
    biasTD500 = mean(modobsTD500(:,1)) - mean(modobsTD500(:,2));
    biasTD700 = mean(modobsTD700(:,1)) - mean(modobsTD700(:,2));
    biasTD850 = mean(modobsTD850(:,1)) - mean(modobsTD850(:,2));
    biasTD925 = mean(modobsTD925(:,1)) - mean(modobsTD925(:,2));
    biasTD1000 = mean(modobsTD1000(:,1)) - mean(modobsTD1000(:,2));
    
    biasWIND150 = mean(modobsWIND150(:,1)) - mean(modobsWIND150(:,2));
    biasWIND200 = mean(modobsWIND200(:,1)) - mean(modobsWIND200(:,2));
    biasWIND250 = mean(modobsWIND250(:,1)) - mean(modobsWIND250(:,2));
    biasWIND300 = mean(modobsWIND300(:,1)) - mean(modobsWIND300(:,2));
    biasWIND400 = mean(modobsWIND400(:,1)) - mean(modobsWIND400(:,2));
    biasWIND500 = mean(modobsWIND500(:,1)) - mean(modobsWIND500(:,2));
    biasWIND700 = mean(modobsWIND700(:,1)) - mean(modobsWIND700(:,2));
    biasWIND850 = mean(modobsWIND850(:,1)) - mean(modobsWIND850(:,2));
    biasWIND925 = mean(modobsWIND925(:,1)) - mean(modobsWIND925(:,2));
    biasWIND1000 = mean(modobsWIND1000(:,1)) - mean(modobsWIND1000(:,2));    
    
    
    rmseT150_all=sqrt((modobsT150(:,1)-modobsT150(:,2)).^2); rmseT150=mean(rmseT150_all);
    rmseT200_all=sqrt((modobsT200(:,1)-modobsT200(:,2)).^2); rmseT200=mean(rmseT200_all);
    rmseT250_all=sqrt((modobsT250(:,1)-modobsT250(:,2)).^2); rmseT250=mean(rmseT250_all);
    rmseT300_all=sqrt((modobsT300(:,1)-modobsT300(:,2)).^2); rmseT300=mean(rmseT300_all);
    rmseT400_all=sqrt((modobsT400(:,1)-modobsT400(:,2)).^2); rmseT400=mean(rmseT400_all);
    rmseT500_all=sqrt((modobsT500(:,1)-modobsT500(:,2)).^2); rmseT500=mean(rmseT500_all);
    rmseT700_all=sqrt((modobsT700(:,1)-modobsT700(:,2)).^2); rmseT700=mean(rmseT700_all);
    rmseT850_all=sqrt((modobsT850(:,1)-modobsT850(:,2)).^2); rmseT850=mean(rmseT850_all);
    rmseT925_all=sqrt((modobsT925(:,1)-modobsT925(:,2)).^2); rmseT925=mean(rmseT925_all);
    rmseT1000_all=sqrt((modobsT1000(:,1)-modobsT1000(:,2)).^2); rmseT1000=mean(rmseT1000_all);

    rmseTD150_all=sqrt((modobsTD150(:,1)-modobsTD150(:,2)).^2); rmseTD150=mean(rmseTD150_all);
    rmseTD200_all=sqrt((modobsTD200(:,1)-modobsTD200(:,2)).^2); rmseTD200=mean(rmseTD200_all);
    rmseTD250_all=sqrt((modobsTD250(:,1)-modobsTD250(:,2)).^2); rmseTD250=mean(rmseTD250_all);
    rmseTD300_all=sqrt((modobsTD300(:,1)-modobsTD300(:,2)).^2); rmseTD300=mean(rmseTD300_all);
    rmseTD400_all=sqrt((modobsTD400(:,1)-modobsTD400(:,2)).^2); rmseTD400=mean(rmseTD400_all);
    rmseTD500_all=sqrt((modobsTD500(:,1)-modobsTD500(:,2)).^2); rmseTD500=mean(rmseTD500_all);
    rmseTD700_all=sqrt((modobsTD700(:,1)-modobsTD700(:,2)).^2); rmseTD700=mean(rmseTD700_all);
    rmseTD850_all=sqrt((modobsTD850(:,1)-modobsTD850(:,2)).^2); rmseTD850=mean(rmseTD850_all);
    rmseTD925_all=sqrt((modobsTD925(:,1)-modobsTD925(:,2)).^2); rmseTD925=mean(rmseTD925_all);
    rmseTD1000_all=sqrt((modobsTD1000(:,1)-modobsTD1000(:,2)).^2); rmseTD1000=mean(rmseTD1000_all);
    
    rmseWIND150_all=sqrt((modobsWIND150(:,1)-modobsWIND150(:,2)).^2); rmseWIND150=mean(rmseWIND150_all);
    rmseWIND200_all=sqrt((modobsWIND200(:,1)-modobsWIND200(:,2)).^2); rmseWIND200=mean(rmseWIND200_all);
    rmseWIND250_all=sqrt((modobsWIND250(:,1)-modobsWIND250(:,2)).^2); rmseWIND250=mean(rmseWIND250_all);
    rmseWIND300_all=sqrt((modobsWIND300(:,1)-modobsWIND300(:,2)).^2); rmseWIND300=mean(rmseWIND300_all);
    rmseWIND400_all=sqrt((modobsWIND400(:,1)-modobsWIND400(:,2)).^2); rmseWIND400=mean(rmseWIND400_all);
    rmseWIND500_all=sqrt((modobsWIND500(:,1)-modobsWIND500(:,2)).^2); rmseWIND500=mean(rmseWIND500_all);
    rmseWIND700_all=sqrt((modobsWIND700(:,1)-modobsWIND700(:,2)).^2); rmseWIND700=mean(rmseWIND700_all);
    rmseWIND850_all=sqrt((modobsWIND850(:,1)-modobsWIND850(:,2)).^2); rmseWIND850=mean(rmseWIND850_all);
    rmseWIND925_all=sqrt((modobsWIND925(:,1)-modobsWIND925(:,2)).^2); rmseWIND925=mean(rmseWIND925_all);
    rmseWIND1000_all=sqrt((modobsWIND1000(:,1)-modobsWIND1000(:,2)).^2); rmseWIND1000=mean(rmseWIND1000_all);
%     rmseWIND150=mean(sqrt((modobsWIND150(:,1)-modobsWIND150(:,2)).^2));
%     rmseWIND200=mean(sqrt((modobsWIND200(:,1)-modobsWIND200(:,2)).^2));
%     rmseWIND250=mean(sqrt((modobsWIND250(:,1)-modobsWIND250(:,2)).^2));
%     rmseWIND300=mean(sqrt((modobsWIND300(:,1)-modobsWIND300(:,2)).^2));
%     rmseWIND400=mean(sqrt((modobsWIND400(:,1)-modobsWIND400(:,2)).^2));
%     rmseWIND500=mean(sqrt((modobsWIND500(:,1)-modobsWIND500(:,2)).^2));
%     rmseWIND700=mean(sqrt((modobsWIND700(:,1)-modobsWIND700(:,2)).^2));
%     rmseWIND850=mean(sqrt((modobsWIND850(:,1)-modobsWIND850(:,2)).^2));
%     rmseWIND925=mean(sqrt((modobsWIND925(:,1)-modobsWIND925(:,2)).^2));
%     rmseWIND1000=mean(sqrt((modobsWIND1000(:,1)-modobsWIND1000(:,2)).^2));



    %%% Rank Histogram, Error to Variance, Variance

    %%%%%% Temperature
    
    for a=1:size(T150,1)
        obs_err_rand=normrnd(0,sigma_T150,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(T150(a,:)+obs_err_rand);
        rank=find(sorted==T150(a,num_mems+1));
        rank_sum_T150(rank,t)=rank_sum_T150(rank,t)+1;
        mod_var(a)=var(T150(a,1:num_mems));
%         stat_errtovarT150(a)=( (modobsT150(a,2)-modobsT150(a,1))' * (modobsT150(a,2)-modobsT150(a,1))) / (mod_var(a)+sigma_T150^2);
    end
    obs_var(1:a)=sigma_T150^2;
    tot_varT150=obs_var+mod_var;
    errtovarT150=( (modobsT150(:,2)-modobsT150(:,1))' * (modobsT150(:,2)-modobsT150(:,1))) / (sum(tot_varT150));
    varT150=mean(mod_var);
    clearvars mod_var obs_var a         
    
    for a=1:size(T200,1)
        obs_err_rand=normrnd(0,sigma_T200,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(T200(a,:)+obs_err_rand);
        rank=find(sorted==T200(a,num_mems+1));
        rank_sum_T200(rank,t)=rank_sum_T200(rank,t)+1;
        mod_var(a)=var(T200(a,1:num_mems));
    end
    obs_var(1:a)=sigma_T200^2;
    tot_varT200=obs_var+mod_var;
    errtovarT200=( (modobsT200(:,2)-modobsT200(:,1))' * (modobsT200(:,2)-modobsT200(:,1))) / (sum(tot_varT200));
    varT200=mean(mod_var);
    clearvars mod_var obs_var a     
    
    for a=1:size(T250,1)
        obs_err_rand=normrnd(0,sigma_T250,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(T250(a,:)+obs_err_rand);
        rank=find(sorted==T250(a,num_mems+1));
        rank_sum_T250(rank,t)=rank_sum_T250(rank,t)+1;
        mod_var(a)=var(T250(a,1:num_mems));
    end
    obs_var(1:a)=sigma_T250^2;
    tot_varT250=obs_var+mod_var;
    errtovarT250=( (modobsT250(:,2)-modobsT250(:,1))' * (modobsT250(:,2)-modobsT250(:,1))) / (sum(tot_varT250));
    varT250=mean(mod_var);
    clearvars mod_var obs_var a    
    
    for a=1:size(T300,1)
        obs_err_rand=normrnd(0,sigma_T300,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(T300(a,:)+obs_err_rand);
        rank=find(sorted==T300(a,num_mems+1));
        rank_sum_T300(rank,t)=rank_sum_T300(rank,t)+1;
        mod_var(a)=var(T300(a,1:num_mems));
    end
    obs_var(1:a)=sigma_T300^2;
    tot_varT300=obs_var+mod_var;
    errtovarT300=( (modobsT300(:,2)-modobsT300(:,1))' * (modobsT300(:,2)-modobsT300(:,1))) / (sum(tot_varT300));
    varT300=mean(mod_var);
    clearvars mod_var obs_var a
    
    for a=1:size(T400,1)
        obs_err_rand=normrnd(0,sigma_T400,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(T400(a,:)+obs_err_rand);
        rank=find(sorted==T400(a,num_mems+1));
        rank_sum_T400(rank,t)=rank_sum_T400(rank,t)+1;
        mod_var(a)=var(T400(a,1:num_mems));
    end
    obs_var(1:a)=sigma_T400^2;
    tot_varT400=obs_var+mod_var;
    errtovarT400=( (modobsT400(:,2)-modobsT400(:,1))' * (modobsT400(:,2)-modobsT400(:,1))) / (sum(tot_varT400));
    varT400=mean(mod_var);
    clearvars mod_var obs_var a 
    
    for a=1:size(T500,1)
        obs_err_rand=normrnd(0,sigma_T500,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(T500(a,:)+obs_err_rand);
        rank=find(sorted==T500(a,num_mems+1));
        rank_sum_T500(rank,t)=rank_sum_T500(rank,t)+1;
        mod_var(a)=var(T500(a,1:num_mems));
    end
    obs_var(1:a)=sigma_T500^2;
    tot_varT500=obs_var+mod_var;
    errtovarT500=( (modobsT500(:,2)-modobsT500(:,1))' * (modobsT500(:,2)-modobsT500(:,1))) / (sum(tot_varT500));
    varT500=mean(mod_var);
    clearvars mod_var obs_var a 
    
    for a=1:size(T700,1)
        obs_err_rand=normrnd(0,sigma_T700,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(T700(a,:)+obs_err_rand);
        rank=find(sorted==T700(a,num_mems+1));
        rank_sum_T700(rank,t)=rank_sum_T700(rank,t)+1;
        mod_var(a)=var(T700(a,1:num_mems));
    end
    obs_var(1:a)=sigma_T700^2;
    tot_varT700=obs_var+mod_var;
    errtovarT700=( (modobsT700(:,2)-modobsT700(:,1))' * (modobsT700(:,2)-modobsT700(:,1))) / (sum(tot_varT700));
    varT700=mean(mod_var);
    clearvars mod_var obs_var a 
    
    for a=1:size(T850,1)
        obs_err_rand=normrnd(0,sigma_T850,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(T850(a,:)+obs_err_rand);
        rank=find(sorted==T850(a,num_mems+1));
        rank_sum_T850(rank,t)=rank_sum_T850(rank,t)+1;
        mod_var(a)=var(T850(a,1:num_mems));
    end
    obs_var(1:a)=sigma_T850^2;
    tot_varT850=obs_var+mod_var;
    errtovarT850=( (modobsT850(:,2)-modobsT850(:,1))' * (modobsT850(:,2)-modobsT850(:,1))) / (sum(tot_varT850));
    varT850=mean(mod_var);
    clearvars mod_var obs_var a 
    
    for a=1:size(T925,1)
        obs_err_rand=normrnd(0,sigma_T925,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(T925(a,:)+obs_err_rand);
        rank=find(sorted==T925(a,num_mems+1));
        rank_sum_T925(rank,t)=rank_sum_T925(rank,t)+1;
        mod_var(a)=var(T925(a,1:num_mems));
    end
    obs_var(1:a)=sigma_T925^2;
    tot_varT925=obs_var+mod_var;
    errtovarT925=( (modobsT925(:,2)-modobsT925(:,1))' * (modobsT925(:,2)-modobsT925(:,1))) / (sum(tot_varT925));
    varT925=mean(mod_var);
    clearvars mod_var obs_var a 

    for a=1:size(T1000,1)
        obs_err_rand=normrnd(0,sigma_T1000,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(T1000(a,:)+obs_err_rand);
        rank=find(sorted==T1000(a,num_mems+1));
        rank_sum_T1000(rank,t)=rank_sum_T1000(rank,t)+1;
        mod_var(a)=var(T1000(a,1:num_mems));
    end
    obs_var(1:a)=sigma_T1000^2;
    tot_varT1000=obs_var+mod_var;
    errtovarT1000=( (modobsT1000(:,2)-modobsT1000(:,1))' * (modobsT1000(:,2)-modobsT1000(:,1))) / (sum(tot_varT1000));
    varT1000=mean(mod_var);
    clearvars mod_var obs_var a     

    %%%% WIND
    for a=1:size(WIND150,1)
        obs_err_rand=normrnd(0,sigma_WIND150,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(WIND150(a,:)+obs_err_rand);
        rank=find(sorted==WIND150(a,num_mems+1));
        rank_sum_WIND150(rank,t)=rank_sum_WIND150(rank,t)+1;
        mod_var(a)=var(WIND150(a,1:num_mems));
    end
    obs_var(1:a)=sigma_WIND150^2;
    tot_varWIND150=obs_var+mod_var;
    errtovarWIND150=( (modobsWIND150(:,2)-modobsWIND150(:,1))' * (modobsWIND150(:,2)-modobsWIND150(:,1))) / (sum(tot_varWIND150));
    varWIND150=mean(mod_var);
    clearvars mod_var obs_var a    
    
     for a=1:size(WIND200,1)
        obs_err_rand=normrnd(0,sigma_WIND200,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(WIND200(a,:)+obs_err_rand);
        rank=find(sorted==WIND200(a,num_mems+1));
        rank_sum_WIND200(rank,t)=rank_sum_WIND200(rank,t)+1;
        mod_var(a)=var(WIND200(a,1:num_mems));
    end
    obs_var(1:a)=sigma_WIND200^2;
    tot_varWIND200=obs_var+mod_var;
    errtovarWIND200=( (modobsWIND200(:,2)-modobsWIND200(:,1))' * (modobsWIND200(:,2)-modobsWIND200(:,1))) / (sum(tot_varWIND200));
    varWIND200=mean(mod_var);
    clearvars mod_var obs_var a   
    
    for a=1:size(WIND250,1)
        obs_err_rand=normrnd(0,sigma_WIND250,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(WIND250(a,:)+obs_err_rand);
        rank=find(sorted==WIND250(a,num_mems+1));
        rank_sum_WIND250(rank,t)=rank_sum_WIND250(rank,t)+1;
        mod_var(a)=var(WIND250(a,1:num_mems));
    end
    obs_var(1:a)=sigma_WIND250^2;
    tot_varWIND250=obs_var+mod_var;
    errtovarWIND250=( (modobsWIND250(:,2)-modobsWIND250(:,1))' * (modobsWIND250(:,2)-modobsWIND250(:,1))) / (sum(tot_varWIND250));
    varWIND250=mean(mod_var);
    clearvars mod_var obs_var a    
    
    for a=1:size(WIND300,1)
        obs_err_rand=normrnd(0,sigma_WIND300,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(WIND300(a,:)+obs_err_rand);
        rank=find(sorted==WIND300(a,num_mems+1));
        rank_sum_WIND300(rank,t)=rank_sum_WIND300(rank,t)+1;
        mod_var(a)=var(WIND300(a,1:num_mems));
    end
    obs_var(1:a)=sigma_WIND300^2;
    tot_varWIND300=obs_var+mod_var;
    errtovarWIND300=( (modobsWIND300(:,2)-modobsWIND300(:,1))' * (modobsWIND300(:,2)-modobsWIND300(:,1))) / (sum(tot_varWIND300));
    varWIND300=mean(mod_var);
    clearvars mod_var obs_var a

    for a=1:size(WIND400,1)
        obs_err_rand=normrnd(0,sigma_WIND400,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(WIND400(a,:)+obs_err_rand);
        rank=find(sorted==WIND400(a,num_mems+1));
        rank_sum_WIND400(rank,t)=rank_sum_WIND400(rank,t)+1;
        mod_var(a)=var(WIND400(a,1:num_mems));
    end
    obs_var(1:a)=sigma_WIND400^2;
    tot_varWIND400=obs_var+mod_var;
    errtovarWIND400=( (modobsWIND400(:,2)-modobsWIND400(:,1))' * (modobsWIND400(:,2)-modobsWIND400(:,1))) / (sum(tot_varWIND400));
    varWIND400=mean(mod_var);
    clearvars mod_var obs_var a
    
    for a=1:size(WIND500,1)
        obs_err_rand=normrnd(0,sigma_WIND500,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(WIND500(a,:)+obs_err_rand);
        rank=find(sorted==WIND500(a,num_mems+1));
        rank_sum_WIND500(rank,t)=rank_sum_WIND500(rank,t)+1;
        mod_var(a)=var(WIND500(a,1:num_mems));
    end
    obs_var(1:a)=sigma_WIND500^2;
    tot_varWIND500=obs_var+mod_var;
    errtovarWIND500=( (modobsWIND500(:,2)-modobsWIND500(:,1))' * (modobsWIND500(:,2)-modobsWIND500(:,1))) / (sum(tot_varWIND500));
    varWIND500=mean(mod_var);
    clearvars mod_var obs_var a    

    for a=1:size(WIND700,1)
        obs_err_rand=normrnd(0,sigma_WIND700,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(WIND700(a,:)+obs_err_rand);
        rank=find(sorted==WIND700(a,num_mems+1));
        rank_sum_WIND700(rank,t)=rank_sum_WIND700(rank,t)+1;
        mod_var(a)=var(WIND700(a,1:num_mems));
    end
    obs_var(1:a)=sigma_WIND700^2;
    tot_varWIND700=obs_var+mod_var;
    errtovarWIND700=( (modobsWIND700(:,2)-modobsWIND700(:,1))' * (modobsWIND700(:,2)-modobsWIND700(:,1))) / (sum(tot_varWIND700));
    varWIND700=mean(mod_var);
    clearvars mod_var obs_var a
    
    for a=1:size(WIND850,1)
        obs_err_rand=normrnd(0,sigma_WIND850,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(WIND850(a,:)+obs_err_rand);
        rank=find(sorted==WIND850(a,num_mems+1));
        rank_sum_WIND850(rank,t)=rank_sum_WIND850(rank,t)+1;
        mod_var(a)=var(WIND850(a,1:num_mems));
    end
    obs_var(1:a)=sigma_WIND850^2;
    tot_varWIND850=obs_var+mod_var;
    errtovarWIND850=( (modobsWIND850(:,2)-modobsWIND850(:,1))' * (modobsWIND850(:,2)-modobsWIND850(:,1))) / (sum(tot_varWIND850));
    varWIND850=mean(mod_var);
    clearvars mod_var obs_var a    
    

    for a=1:size(WIND925,1)
        obs_err_rand=normrnd(0,sigma_WIND925,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(WIND925(a,:)+obs_err_rand);
        rank=find(sorted==WIND925(a,num_mems+1));
        rank_sum_WIND925(rank,t)=rank_sum_WIND925(rank,t)+1;
        mod_var(a)=var(WIND925(a,1:num_mems));
    end
    obs_var(1:a)=sigma_WIND925^2;
    tot_varWIND925=obs_var+mod_var;
    errtovarWIND925=( (modobsWIND925(:,2)-modobsWIND925(:,1))' * (modobsWIND925(:,2)-modobsWIND925(:,1))) / (sum(tot_varWIND925));
    varWIND925=mean(mod_var);
    clearvars mod_var obs_var a
    
    for a=1:size(WIND1000,1)
        obs_err_rand=normrnd(0,sigma_WIND1000,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(WIND1000(a,:)+obs_err_rand);
        rank=find(sorted==WIND1000(a,num_mems+1));
        rank_sum_WIND1000(rank,t)=rank_sum_WIND1000(rank,t)+1;
        mod_var(a)=var(WIND1000(a,1:num_mems));
    end
    obs_var(1:a)=sigma_WIND1000^2;
    tot_varWIND1000=obs_var+mod_var;
    errtovarWIND1000=( (modobsWIND1000(:,2)-modobsWIND1000(:,1))' * (modobsWIND1000(:,2)-modobsWIND1000(:,1))) / (sum(tot_varWIND1000));
    varWIND1000=mean(mod_var);
    clearvars mod_var obs_var a    

    
    
    %%%%%% Dewpoint Temperature
   
    for a=1:size(TD150,1)
        obs_err_rand=normrnd(0,sigma_TD150,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(TD150(a,:)+obs_err_rand);
        rank=find(sorted==TD150(a,num_mems+1));
        rank_sum_TD150(rank,t)=rank_sum_TD150(rank,t)+1;
        mod_var(a)=var(TD150(a,1:num_mems));
    end
    obs_var(1:a)=sigma_TD150^2;
    tot_varTD150=obs_var+mod_var;
    errtovarTD150=( (modobsTD150(:,2)-modobsTD150(:,1))' * (modobsTD150(:,2)-modobsTD150(:,1))) / (sum(tot_varTD150));
    varTD150=mean(mod_var);
    clearvars mod_var obs_var a         
    
    for a=1:size(TD200,1)
        obs_err_rand=normrnd(0,sigma_TD200,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(TD200(a,:)+obs_err_rand);
        rank=find(sorted==TD200(a,num_mems+1));
        rank_sum_TD200(rank,t)=rank_sum_TD200(rank,t)+1;
        mod_var(a)=var(TD200(a,1:num_mems));
    end
    obs_var(1:a)=sigma_TD200^2;
    tot_varTD200=obs_var+mod_var;
    errtovarTD200=( (modobsTD200(:,2)-modobsTD200(:,1))' * (modobsTD200(:,2)-modobsTD200(:,1))) / (sum(tot_varTD200));
    varTD200=mean(mod_var);
    clearvars mod_var obs_var a     
    
    for a=1:size(TD250,1)
        obs_err_rand=normrnd(0,sigma_TD250,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(TD250(a,:)+obs_err_rand);
        rank=find(sorted==TD250(a,num_mems+1));
        rank_sum_TD250(rank,t)=rank_sum_TD250(rank,t)+1;
        mod_var(a)=var(TD250(a,1:num_mems));
    end
    obs_var(1:a)=sigma_TD250^2;
    tot_varTD250=obs_var+mod_var;
    errtovarTD250=( (modobsTD250(:,2)-modobsTD250(:,1))' * (modobsTD250(:,2)-modobsTD250(:,1))) / (sum(tot_varTD250));
    varTD250=mean(mod_var);
    clearvars mod_var obs_var a    
    
    for a=1:size(TD300,1)
        obs_err_rand=normrnd(0,sigma_TD300,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(TD300(a,:)+obs_err_rand);
        rank=find(sorted==TD300(a,num_mems+1));
        rank_sum_TD300(rank,t)=rank_sum_TD300(rank,t)+1;
        mod_var(a)=var(TD300(a,1:num_mems));
    end
    obs_var(1:a)=sigma_TD300^2;
    tot_varTD300=obs_var+mod_var;
    errtovarTD300=( (modobsTD300(:,2)-modobsTD300(:,1))' * (modobsTD300(:,2)-modobsTD300(:,1))) / (sum(tot_varTD300));
    varTD300=mean(mod_var);
    clearvars mod_var obs_var a
    
    for a=1:size(TD400,1)
        obs_err_rand=normrnd(0,sigma_TD400,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(TD400(a,:)+obs_err_rand);
        rank=find(sorted==TD400(a,num_mems+1));
        rank_sum_TD400(rank,t)=rank_sum_TD400(rank,t)+1;
        mod_var(a)=var(TD400(a,1:num_mems));
    end
    obs_var(1:a)=sigma_TD400^2;
    tot_varTD400=obs_var+mod_var;
    errtovarTD400=( (modobsTD400(:,2)-modobsTD400(:,1))' * (modobsTD400(:,2)-modobsTD400(:,1))) / (sum(tot_varTD400));
    varTD400=mean(mod_var);
    clearvars mod_var obs_var a 
    
    for a=1:size(TD500,1)
        obs_err_rand=normrnd(0,sigma_TD500,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(TD500(a,:)+obs_err_rand);
        rank=find(sorted==TD500(a,num_mems+1));
        rank_sum_TD500(rank,t)=rank_sum_TD500(rank,t)+1;
        mod_var(a)=var(TD500(a,1:num_mems));
    end
    obs_var(1:a)=sigma_TD500^2;
    tot_varTD500=obs_var+mod_var;
    errtovarTD500=( (modobsTD500(:,2)-modobsTD500(:,1))' * (modobsTD500(:,2)-modobsTD500(:,1))) / (sum(tot_varTD500));
    varTD500=mean(mod_var);
    clearvars mod_var obs_var a 
    
    for a=1:size(TD700,1)
        obs_err_rand=normrnd(0,sigma_TD700,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(TD700(a,:)+obs_err_rand);
        rank=find(sorted==TD700(a,num_mems+1));
        rank_sum_TD700(rank,t)=rank_sum_TD700(rank,t)+1;
        mod_var(a)=var(TD700(a,1:num_mems));
    end
    obs_var(1:a)=sigma_TD700^2;
    tot_varTD700=obs_var+mod_var;
    errtovarTD700=( (modobsTD700(:,2)-modobsTD700(:,1))' * (modobsTD700(:,2)-modobsTD700(:,1))) / (sum(tot_varTD700));
    varTD700=mean(mod_var);
    clearvars mod_var obs_var a 
    
    for a=1:size(TD850,1)
        obs_err_rand=normrnd(0,sigma_TD850,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(TD850(a,:)+obs_err_rand);
        rank=find(sorted==TD850(a,num_mems+1));
        rank_sum_TD850(rank,t)=rank_sum_TD850(rank,t)+1;
        mod_var(a)=var(TD850(a,1:num_mems));
    end
    obs_var(1:a)=sigma_TD850^2;
    tot_varTD850=obs_var+mod_var;
    errtovarTD850=( (modobsTD850(:,2)-modobsTD850(:,1))' * (modobsTD850(:,2)-modobsTD850(:,1))) / (sum(tot_varTD850));
    varTD850=mean(mod_var);
    clearvars mod_var obs_var a 
    
    for a=1:size(TD925,1)
        obs_err_rand=normrnd(0,sigma_TD925,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(TD925(a,:)+obs_err_rand);
        rank=find(sorted==TD925(a,num_mems+1));
        rank_sum_TD925(rank,t)=rank_sum_TD925(rank,t)+1;
        mod_var(a)=var(TD925(a,1:num_mems));
    end
    obs_var(1:a)=sigma_TD925^2;
    tot_varTD925=obs_var+mod_var;
    errtovarTD925=( (modobsTD925(:,2)-modobsTD925(:,1))' * (modobsTD925(:,2)-modobsTD925(:,1))) / (sum(tot_varTD925));
    varTD925=mean(mod_var);
    clearvars mod_var obs_var a 

    for a=1:size(TD1000,1)
        obs_err_rand=normrnd(0,sigma_TD1000,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(TD1000(a,:)+obs_err_rand);
        rank=find(sorted==TD1000(a,num_mems+1));
        rank_sum_TD1000(rank,t)=rank_sum_TD1000(rank,t)+1;
        mod_var(a)=var(TD1000(a,1:num_mems));
    end
    obs_var(1:a)=sigma_TD1000^2;
    tot_varTD1000=obs_var+mod_var;
    errtovarTD1000=( (modobsTD1000(:,2)-modobsTD1000(:,1))' * (modobsTD1000(:,2)-modobsTD1000(:,1))) / (sum(tot_varTD1000));
    varTD1000=mean(mod_var);
    clearvars mod_var obs_var a     
        
   
   
       
    
    
% end

% cd('..')
% end

%Sum Rank Hists through all time periods

if strcmp(type,'fix')==1
%%% FIXED
% fixed_rank_sum_GPH300_accum=sum(rank_sum_GPH300,2);
% fixed_rank_sum_GPH400_accum=sum(rank_sum_GPH400,2);
% fixed_rank_sum_GPH500_accum=sum(rank_sum_GPH500,2);
% fixed_rank_sum_GPH700_accum=sum(rank_sum_GPH700,2);
% fixed_rank_sum_GPH850_accum=sum(rank_sum_GPH850,2);
% fixed_rank_sum_GPH925_accum=sum(rank_sum_GPH925,2);

fixed_rank_sum_T300_accum=sum(rank_sum_T300,2);
fixed_rank_sum_T400_accum=sum(rank_sum_T400,2);
fixed_rank_sum_T500_accum=sum(rank_sum_T500,2);
fixed_rank_sum_T700_accum=sum(rank_sum_T700,2);
fixed_rank_sum_T850_accum=sum(rank_sum_T850,2);
fixed_rank_sum_T925_accum=sum(rank_sum_T925,2);

fixed_rank_sum_WIND300_accum=sum(rank_sum_WIND300,2);
fixed_rank_sum_WIND400_accum=sum(rank_sum_WIND400,2);
fixed_rank_sum_WIND500_accum=sum(rank_sum_WIND500,2);
fixed_rank_sum_WIND700_accum=sum(rank_sum_WIND700,2);
fixed_rank_sum_WIND850_accum=sum(rank_sum_WIND850,2);
fixed_rank_sum_WIND925_accum=sum(rank_sum_WIND925,2);

elseif strcmp(type,'var')==1
%%% VARIED
% varied_rank_sum_GPH300_accum=sum(rank_sum_GPH300,2);
% varied_rank_sum_GPH400_accum=sum(rank_sum_GPH400,2);
% varied_rank_sum_GPH500_accum=sum(rank_sum_GPH500,2);
% varied_rank_sum_GPH700_accum=sum(rank_sum_GPH700,2);
% varied_rank_sum_GPH850_accum=sum(rank_sum_GPH850,2);
% varied_rank_sum_GPH925_accum=sum(rank_sum_GPH925,2);

varied_rank_sum_T300_accum=sum(rank_sum_T300,2);
varied_rank_sum_T400_accum=sum(rank_sum_T400,2);
varied_rank_sum_T500_accum=sum(rank_sum_T500,2);
varied_rank_sum_T700_accum=sum(rank_sum_T700,2);
varied_rank_sum_T850_accum=sum(rank_sum_T850,2);
varied_rank_sum_T925_accum=sum(rank_sum_T925,2);

varied_rank_sum_WIND300_accum=sum(rank_sum_WIND300,2);
varied_rank_sum_WIND400_accum=sum(rank_sum_WIND400,2);
varied_rank_sum_WIND500_accum=sum(rank_sum_WIND500,2);
varied_rank_sum_WIND700_accum=sum(rank_sum_WIND700,2);
varied_rank_sum_WIND850_accum=sum(rank_sum_WIND850,2);
varied_rank_sum_WIND925_accum=sum(rank_sum_WIND925,2);
end
% end
% %%% DIFF (VAR-FIX)
% diff_rank_sum_GPH300_accum=varied_rank_sum_GPH300_accum-fixed_rank_sum_GPH300_accum;
% diff_rank_sum_GPH400_accum=varied_rank_sum_GPH400_accum-fixed_rank_sum_GPH400_accum;
% diff_rank_sum_GPH500_accum=varied_rank_sum_GPH500_accum-fixed_rank_sum_GPH500_accum;
% diff_rank_sum_GPH700_accum=varied_rank_sum_GPH700_accum-fixed_rank_sum_GPH700_accum;
% diff_rank_sum_GPH850_accum=varied_rank_sum_GPH850_accum-fixed_rank_sum_GPH850_accum;
% diff_rank_sum_GPH925_accum=varied_rank_sum_GPH925_accum-fixed_rank_sum_GPH925_accum;
% 
% diff_rank_sum_T300_accum=varied_rank_sum_T300_accum-fixed_rank_sum_T300_accum;
% diff_rank_sum_T400_accum=varied_rank_sum_T400_accum-fixed_rank_sum_T400_accum;
% diff_rank_sum_T500_accum=varied_rank_sum_T500_accum-fixed_rank_sum_T500_accum;
% diff_rank_sum_T700_accum=varied_rank_sum_T700_accum-fixed_rank_sum_T700_accum;
% diff_rank_sum_T850_accum=varied_rank_sum_T850_accum-fixed_rank_sum_T850_accum;
% diff_rank_sum_T925_accum=varied_rank_sum_T925_accum-fixed_rank_sum_T925_accum;
% 
% diff_rank_sum_WIND300_accum=varied_rank_sum_WIND300_accum-fixed_rank_sum_WIND300_accum;
% diff_rank_sum_WIND400_accum=varied_rank_sum_WIND400_accum-fixed_rank_sum_WIND400_accum;
% diff_rank_sum_WIND500_accum=varied_rank_sum_WIND500_accum-fixed_rank_sum_WIND500_accum;
% diff_rank_sum_WIND700_accum=varied_rank_sum_WIND700_accum-fixed_rank_sum_WIND700_accum;
% diff_rank_sum_WIND850_accum=varied_rank_sum_WIND850_accum-fixed_rank_sum_WIND850_accum;
% diff_rank_sum_WIND925_accum=varied_rank_sum_WIND925_accum-fixed_rank_sum_WIND925_accum;


% figure
% bar(rank_sum_GPH300_accum)
% xlim([1 51])
%


%%%%%%%%%%%%%%%%%%%%%%% RANK PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %T
% figure
% subplot(6,1,1); bar(eval([type '_rank_sum_T300_accum'])); xlim([1 51]); %ylim([0 8])
% title('300 hPa temperature')
% ylabel('Count')
% subplot(6,1,2); bar(eval([type '_rank_sum_T400_accum'])); xlim([1 51]); %ylim([0 8])
% title('400 hPa temperature')
% ylabel('Count')
% subplot(6,1,3); bar(eval([type '_rank_sum_T500_accum'])); xlim([1 51]); %ylim([0 8])
% title('500 hPa temperature')
% ylabel('Count')
% subplot(6,1,4); bar(eval([type '_rank_sum_T700_accum'])); xlim([1 51]); %ylim([0 8])
% title('700 hPa temperature')
% ylabel('Count')
% subplot(6,1,5); bar(eval([type '_rank_sum_T850_accum'])); xlim([1 51]); %ylim([0 8])
% title('850 hPa temperature')
% ylabel('Count')
% subplot(6,1,6); bar(eval([type '_rank_sum_T925_accum'])); xlim([1 51]); %ylim([0 8])
% title('925 hPa temperature')
% ylabel('Count')
% xlabel('Observation Rank')
% %WIND
% figure
% subplot(6,1,1); bar(eval([type '_rank_sum_WIND300_accum'])); xlim([1 51]); %ylim([0 8])
% title('300 hPa horizontal wind')
% ylabel('Count')
% subplot(6,1,2); bar(eval([type '_rank_sum_WIND400_accum'])); xlim([1 51]); %ylim([0 8])
% title('400 hPa horizontal wind')
% ylabel('Count')
% subplot(6,1,3); bar(eval([type '_rank_sum_WIND500_accum'])); xlim([1 51]); %ylim([0 8])
% title('500 hPa horizontal wind')
% ylabel('Count')
% subplot(6,1,4); bar(eval([type '_rank_sum_WIND700_accum'])); xlim([1 51]); %ylim([0 8])
% title('700 hPa horizontal wind')
% ylabel('Count')
% subplot(6,1,5); bar(eval([type '_rank_sum_WIND850_accum'])); xlim([1 51]); %ylim([0 8])
% title('850 hPa horizontal wind')
% ylabel('Count')
% subplot(6,1,6); bar(eval([type '_rank_sum_WIND925_accum'])); xlim([1 51]); %ylim([0 8])
% title('925 hPa horizontal wind')
% ylabel('Count')
% xlabel('Observation Rank')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% var='T';
% etvheight_T(6)=mean(mean(eval(['errtovar' var '300'])));
% etvheight_T(5)=mean(mean(eval(['errtovar' var '400'])));
% etvheight_T(4)=mean(mean(eval(['errtovar' var '500'])));
% etvheight_T(3)=mean(mean(eval(['errtovar' var '700'])));
% etvheight_T(2)=mean(mean(eval(['errtovar' var '850'])));
% etvheight_T(1)=mean(mean(eval(['errtovar' var '925'])));
% var='WIND';
% etvheight_WIND(6)=mean(mean(eval(['errtovar' var '300'])));
% etvheight_WIND(5)=mean(mean(eval(['errtovar' var '400'])));
% etvheight_WIND(4)=mean(mean(eval(['errtovar' var '500'])));
% etvheight_WIND(3)=mean(mean(eval(['errtovar' var '700'])));
% etvheight_WIND(2)=mean(mean(eval(['errtovar' var '850'])));
% etvheight_WIND(1)=mean(mean(eval(['errtovar' var '925'])));

% height=[925,850,700,500,400,300];
% ideal=[1,1,1,1,1,1];
% figure
% subplot(1,2,1)
% plot(etvheight_WIND,height,'b'); hold on
% set(gca,'ydir','reverse')
% hold on; plot(ideal,height,'k')
% legend('D2D3')
% title('Horizontal Wind (m/s)'); hold on; axis tight; xlim([0 15])
% ylabel('Pressure (hPa)')
% xlabel('Error-to-Variance')
% 
% subplot(1,2,2)
% set(gca,'ydir','reverse')
% plot(etvheight_T,height,'b'); hold on
% set(gca,'ydir','reverse')
% hold on; plot(ideal,height,'k')
% legend('D2D3')
% title('Temperature (K)'); hold on; axis tight; xlim([0 15])
% ylabel('Pressure (hPa)')
% xlabel('Error-to-Variance')
% total_T150 = sum(sum(obs_count_T150));
 

var='T';
varianceheight_T(10)=eval(['var' var '150']);
varianceheight_T(9)=eval(['var' var '200']);
varianceheight_T(8)=eval(['var' var '250']);
varianceheight_T(7)=eval(['var' var '300']);
varianceheight_T(6)=eval(['var' var '400']);
varianceheight_T(5)=eval(['var' var '500']);
varianceheight_T(4)=eval(['var' var '700']);
varianceheight_T(3)=eval(['var' var '850']);
varianceheight_T(2)=eval(['var' var '925']);
varianceheight_T(1)=eval(['var' var '1000']);

rmseheight_T(10)=eval(['rmse' var '150']);
rmseheight_T(9)=eval(['rmse' var '200']);
rmseheight_T(8)=eval(['rmse' var '250']);
rmseheight_T(7)=eval(['rmse' var '300']);
rmseheight_T(6)=eval(['rmse' var '400']);
rmseheight_T(5)=eval(['rmse' var '500']);
rmseheight_T(4)=eval(['rmse' var '700']);
rmseheight_T(3)=eval(['rmse' var '850']);
rmseheight_T(2)=eval(['rmse' var '925']);
rmseheight_T(1)=eval(['rmse' var '1000']);

etvheight_T(10)=eval(['errtovar' var '150']);
etvheight_T(9)=eval(['errtovar' var '200']);
etvheight_T(8)=eval(['errtovar' var '250']);
etvheight_T(7)=eval(['errtovar' var '300']);
etvheight_T(6)=eval(['errtovar' var '400']);
etvheight_T(5)=eval(['errtovar' var '500']);
etvheight_T(4)=eval(['errtovar' var '700']);
etvheight_T(3)=eval(['errtovar' var '850']);
etvheight_T(2)=eval(['errtovar' var '925']);
etvheight_T(1)=eval(['errtovar' var '1000']);

biasheight_T(10)=eval(['bias' var '150']);
biasheight_T(9)=eval(['bias' var '200']);
biasheight_T(8)=eval(['bias' var '250']);
biasheight_T(7)=eval(['bias' var '300']);
biasheight_T(6)=eval(['bias' var '400']);
biasheight_T(5)=eval(['bias' var '500']);
biasheight_T(4)=eval(['bias' var '700']);
biasheight_T(3)=eval(['bias' var '850']);
biasheight_T(2)=eval(['bias' var '925']);
biasheight_T(1)=eval(['bias' var '1000']);





var='TD';
varianceheight_TD(10)=eval(['var' var '150']);
varianceheight_TD(9)=eval(['var' var '200']);
varianceheight_TD(8)=eval(['var' var '250']);
varianceheight_TD(7)=eval(['var' var '300']);
varianceheight_TD(6)=eval(['var' var '400']);
varianceheight_TD(5)=eval(['var' var '500']);
varianceheight_TD(4)=eval(['var' var '700']);
varianceheight_TD(3)=eval(['var' var '850']);
varianceheight_TD(2)=eval(['var' var '925']);
varianceheight_TD(1)=eval(['var' var '1000']);

rmseheight_TD(10)=eval(['rmse' var '150']);
rmseheight_TD(9)=eval(['rmse' var '200']);
rmseheight_TD(8)=eval(['rmse' var '250']);
rmseheight_TD(7)=eval(['rmse' var '300']);
rmseheight_TD(6)=eval(['rmse' var '400']);
rmseheight_TD(5)=eval(['rmse' var '500']);
rmseheight_TD(4)=eval(['rmse' var '700']);
rmseheight_TD(3)=eval(['rmse' var '850']);
rmseheight_TD(2)=eval(['rmse' var '925']);
rmseheight_TD(1)=eval(['rmse' var '1000']);

etvheight_TD(10)=eval(['errtovar' var '150']);
etvheight_TD(9)=eval(['errtovar' var '200']);
etvheight_TD(8)=eval(['errtovar' var '250']);
etvheight_TD(7)=eval(['errtovar' var '300']);
etvheight_TD(6)=eval(['errtovar' var '400']);
etvheight_TD(5)=eval(['errtovar' var '500']);
etvheight_TD(4)=eval(['errtovar' var '700']);
etvheight_TD(3)=eval(['errtovar' var '850']);
etvheight_TD(2)=eval(['errtovar' var '925']);
etvheight_TD(1)=eval(['errtovar' var '1000']);

biasheight_TD(10)=eval(['bias' var '150']);
biasheight_TD(9)=eval(['bias' var '200']);
biasheight_TD(8)=eval(['bias' var '250']);
biasheight_TD(7)=eval(['bias' var '300']);
biasheight_TD(6)=eval(['bias' var '400']);
biasheight_TD(5)=eval(['bias' var '500']);
biasheight_TD(4)=eval(['bias' var '700']);
biasheight_TD(3)=eval(['bias' var '850']);
biasheight_TD(2)=eval(['bias' var '925']);
biasheight_TD(1)=eval(['bias' var '1000']);





var='WIND';
varianceheight_WIND(10)=eval(['var' var '150']);
varianceheight_WIND(9)=eval(['var' var '200']);
varianceheight_WIND(8)=eval(['var' var '250']);
varianceheight_WIND(7)=eval(['var' var '300']);
varianceheight_WIND(6)=eval(['var' var '400']);
varianceheight_WIND(5)=eval(['var' var '500']);
varianceheight_WIND(4)=eval(['var' var '700']);
varianceheight_WIND(3)=eval(['var' var '850']);
varianceheight_WIND(2)=eval(['var' var '925']);
varianceheight_WIND(1)=eval(['var' var '1000']);

rmseheight_WIND(10)=eval(['rmse' var '150']);
rmseheight_WIND(9)=eval(['rmse' var '200']);
rmseheight_WIND(8)=eval(['rmse' var '250']);
rmseheight_WIND(7)=eval(['rmse' var '300']);
rmseheight_WIND(6)=eval(['rmse' var '400']);
rmseheight_WIND(5)=eval(['rmse' var '500']);
rmseheight_WIND(4)=eval(['rmse' var '700']);
rmseheight_WIND(3)=eval(['rmse' var '850']);
rmseheight_WIND(2)=eval(['rmse' var '925']);
rmseheight_WIND(1)=eval(['rmse' var '1000']);

etvheight_WIND(10)=eval(['errtovar' var '150']);
etvheight_WIND(9)=eval(['errtovar' var '200']);
etvheight_WIND(8)=eval(['errtovar' var '250']);
etvheight_WIND(7)=eval(['errtovar' var '300']);
etvheight_WIND(6)=eval(['errtovar' var '400']);
etvheight_WIND(5)=eval(['errtovar' var '500']);
etvheight_WIND(4)=eval(['errtovar' var '700']);
etvheight_WIND(3)=eval(['errtovar' var '850']);
etvheight_WIND(2)=eval(['errtovar' var '925']);
etvheight_WIND(1)=eval(['errtovar' var '1000']);

biasheight_WIND(10)=eval(['bias' var '150']);
biasheight_WIND(9)=eval(['bias' var '200']);
biasheight_WIND(8)=eval(['bias' var '250']);
biasheight_WIND(7)=eval(['bias' var '300']);
biasheight_WIND(6)=eval(['bias' var '400']);
biasheight_WIND(5)=eval(['bias' var '500']);
biasheight_WIND(4)=eval(['bias' var '700']);
biasheight_WIND(3)=eval(['bias' var '850']);
biasheight_WIND(2)=eval(['bias' var '925']);
biasheight_WIND(1)=eval(['bias' var '1000']);




disp(['domain is: ' domain ' and physics type is: ' type ])

% clearvars -except varianceheight* rmseheight* etvheight*
%%%%%%%%%%%%%% E/V PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(start)

%%% BOOTSTRAP RESAMPLING %%%
if strcmp(type,'fix') == 1 && strcmp(typed1,'fix') == 1  

resample = 1000;

%%% resample the error and variance for each variable, level
rs_etvT150 = zeros(1,resample);
samplelength = length(rmseT150_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errT150 = rmseT150_all(resample_inds).^2;   rs_varT150 = tot_varT150(resample_inds);

rs_etvT200 = zeros(1,resample);
samplelength = length(rmseT200_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errT200 = rmseT200_all(resample_inds).^2;   rs_varT200 = tot_varT200(resample_inds);

rs_etvT250 = zeros(1,resample);
samplelength = length(rmseT250_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errT250 = rmseT250_all(resample_inds).^2;   rs_varT250 = tot_varT250(resample_inds);

rs_etvT300 = zeros(1,resample);
samplelength = length(rmseT300_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errT300 = rmseT300_all(resample_inds).^2;   rs_varT300 = tot_varT300(resample_inds);

rs_etvT400 = zeros(1,resample);
samplelength = length(rmseT400_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errT400 = rmseT400_all(resample_inds).^2;   rs_varT400 = tot_varT400(resample_inds);

rs_etvT500 = zeros(1,resample);
samplelength = length(rmseT500_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errT500 = rmseT500_all(resample_inds).^2;   rs_varT500 = tot_varT500(resample_inds);

rs_etvT700 = zeros(1,resample);
samplelength = length(rmseT700_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errT700 = rmseT700_all(resample_inds).^2;   rs_varT700 = tot_varT700(resample_inds);

rs_etvT850 = zeros(1,resample);
samplelength = length(rmseT850_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errT850 = rmseT850_all(resample_inds).^2;   rs_varT850 = tot_varT850(resample_inds);

rs_etvT925 = zeros(1,resample);
samplelength = length(rmseT925_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errT925 = rmseT925_all(resample_inds).^2;   rs_varT925 = tot_varT925(resample_inds);

rs_etvT1000 = zeros(1,resample);
samplelength = length(rmseT1000_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errT1000 = rmseT1000_all(resample_inds).^2;   rs_varT1000 = tot_varT1000(resample_inds);



rs_etvTD150 = zeros(1,resample);
samplelength = length(rmseTD150_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errTD150 = rmseTD150_all(resample_inds).^2;   rs_varTD150 = tot_varTD150(resample_inds);

rs_etvTD200 = zeros(1,resample);
samplelength = length(rmseTD200_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errTD200 = rmseTD200_all(resample_inds).^2;   rs_varTD200 = tot_varTD200(resample_inds);

rs_etvTD250 = zeros(1,resample);
samplelength = length(rmseTD250_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errTD250 = rmseTD250_all(resample_inds).^2;   rs_varTD250 = tot_varTD250(resample_inds);

rs_etvTD300 = zeros(1,resample);
samplelength = length(rmseTD300_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errTD300 = rmseTD300_all(resample_inds).^2;   rs_varTD300 = tot_varTD300(resample_inds);

rs_etvTD400 = zeros(1,resample);
samplelength = length(rmseTD400_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errTD400 = rmseTD400_all(resample_inds).^2;   rs_varTD400 = tot_varTD400(resample_inds);

rs_etvTD500 = zeros(1,resample);
samplelength = length(rmseTD500_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errTD500 = rmseTD500_all(resample_inds).^2;   rs_varTD500 = tot_varTD500(resample_inds);

rs_etvTD700 = zeros(1,resample);
samplelength = length(rmseTD700_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errTD700 = rmseTD700_all(resample_inds).^2;   rs_varTD700 = tot_varTD700(resample_inds);

rs_etvTD850 = zeros(1,resample);
samplelength = length(rmseTD850_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errTD850 = rmseTD850_all(resample_inds).^2;   rs_varTD850 = tot_varTD850(resample_inds);

rs_etvTD925 = zeros(1,resample);
samplelength = length(rmseTD925_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errTD925 = rmseTD925_all(resample_inds).^2;   rs_varTD925 = tot_varTD925(resample_inds);

rs_etvTD1000 = zeros(1,resample);
samplelength = length(rmseTD1000_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errTD1000 = rmseTD1000_all(resample_inds).^2;   rs_varTD1000 = tot_varTD1000(resample_inds);


rs_etvWIND150 = zeros(1,resample);
samplelength = length(rmseWIND150_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errWIND150 = rmseWIND150_all(resample_inds).^2;   rs_varWIND150 = tot_varWIND150(resample_inds);

rs_etvWIND200 = zeros(1,resample);
samplelength = length(rmseWIND200_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errWIND200 = rmseWIND200_all(resample_inds).^2;   rs_varWIND200 = tot_varWIND200(resample_inds);

rs_etvWIND250 = zeros(1,resample);
samplelength = length(rmseWIND250_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errWIND250 = rmseWIND250_all(resample_inds).^2;   rs_varWIND250 = tot_varWIND250(resample_inds);

rs_etvWIND300 = zeros(1,resample);
samplelength = length(rmseWIND300_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errWIND300 = rmseWIND300_all(resample_inds).^2;   rs_varWIND300 = tot_varWIND300(resample_inds);

rs_etvWIND400 = zeros(1,resample);
samplelength = length(rmseWIND400_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errWIND400 = rmseWIND400_all(resample_inds).^2;   rs_varWIND400 = tot_varWIND400(resample_inds);

rs_etvWIND500 = zeros(1,resample);
samplelength = length(rmseWIND500_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errWIND500 = rmseWIND500_all(resample_inds).^2;   rs_varWIND500 = tot_varWIND500(resample_inds);

rs_etvWIND700 = zeros(1,resample);
samplelength = length(rmseWIND700_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errWIND700 = rmseWIND700_all(resample_inds).^2;   rs_varWIND700 = tot_varWIND700(resample_inds);

rs_etvWIND850 = zeros(1,resample);
samplelength = length(rmseWIND850_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errWIND850 = rmseWIND850_all(resample_inds).^2;   rs_varWIND850 = tot_varWIND850(resample_inds);

rs_etvWIND925 = zeros(1,resample);
samplelength = length(rmseWIND925_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errWIND925 = rmseWIND925_all(resample_inds).^2;   rs_varWIND925 = tot_varWIND925(resample_inds);

rs_etvWIND1000 = zeros(1,resample);
samplelength = length(rmseWIND1000_all);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errWIND1000 = rmseWIND1000_all(resample_inds).^2;   rs_varWIND1000 = tot_varWIND1000(resample_inds);


for i=1:resample;
    % height=[1000,925,850,700,500,400,300,250,200,150];
    rs_etvT150(i) = sum(rs_errT150(:,i)) / sum(rs_varT150(:,i));
    rs_etvT200(i) = sum(rs_errT200(:,i)) / sum(rs_varT200(:,i));
    rs_etvT250(i) = sum(rs_errT250(:,i)) / sum(rs_varT250(:,i));    
    rs_etvT300(i) = sum(rs_errT300(:,i)) / sum(rs_varT300(:,i));
    rs_etvT400(i) = sum(rs_errT400(:,i)) / sum(rs_varT400(:,i));
    rs_etvT500(i) = sum(rs_errT500(:,i)) / sum(rs_varT500(:,i));
    rs_etvT700(i) = sum(rs_errT700(:,i)) / sum(rs_varT700(:,i));
    rs_etvT850(i) = sum(rs_errT850(:,i)) / sum(rs_varT850(:,i));
    rs_etvT925(i) = sum(rs_errT925(:,i)) / sum(rs_varT925(:,i));
    rs_etvT1000(i) = sum(rs_errT1000(:,i)) / sum(rs_varT1000(:,i));
    
    
    rs_etvTD150(i) = sum(rs_errTD150(:,i)) / sum(rs_varTD150(:,i));
    rs_etvTD200(i) = sum(rs_errTD200(:,i)) / sum(rs_varTD200(:,i));
    rs_etvTD250(i) = sum(rs_errTD250(:,i)) / sum(rs_varTD250(:,i));    
    rs_etvTD300(i) = sum(rs_errTD300(:,i)) / sum(rs_varTD300(:,i));
    rs_etvTD400(i) = sum(rs_errTD400(:,i)) / sum(rs_varTD400(:,i));
    rs_etvTD500(i) = sum(rs_errTD500(:,i)) / sum(rs_varTD500(:,i));
    rs_etvTD700(i) = sum(rs_errTD700(:,i)) / sum(rs_varTD700(:,i));
    rs_etvTD850(i) = sum(rs_errTD850(:,i)) / sum(rs_varTD850(:,i));
    rs_etvTD925(i) = sum(rs_errTD925(:,i)) / sum(rs_varTD925(:,i));
    rs_etvTD1000(i) = sum(rs_errTD1000(:,i)) / sum(rs_varTD1000(:,i)); 
    
    
    rs_etvWIND150(i) = sum(rs_errWIND150(:,i)) / sum(rs_varWIND150(:,i));
    rs_etvWIND200(i) = sum(rs_errWIND200(:,i)) / sum(rs_varWIND200(:,i));
    rs_etvWIND250(i) = sum(rs_errWIND250(:,i)) / sum(rs_varWIND250(:,i));    
    rs_etvWIND300(i) = sum(rs_errWIND300(:,i)) / sum(rs_varWIND300(:,i));
    rs_etvWIND400(i) = sum(rs_errWIND400(:,i)) / sum(rs_varWIND400(:,i));
    rs_etvWIND500(i) = sum(rs_errWIND500(:,i)) / sum(rs_varWIND500(:,i));
    rs_etvWIND700(i) = sum(rs_errWIND700(:,i)) / sum(rs_varWIND700(:,i));
    rs_etvWIND850(i) = sum(rs_errWIND850(:,i)) / sum(rs_varWIND850(:,i));
    rs_etvWIND925(i) = sum(rs_errWIND925(:,i)) / sum(rs_varWIND925(:,i));
    rs_etvWIND1000(i) = sum(rs_errWIND1000(:,i)) / sum(rs_varWIND1000(:,i));    
%     
end




boot_etvT150 = sort(rs_etvT150);
boot_etvT200 = sort(rs_etvT200);
boot_etvT250 = sort(rs_etvT250);
boot_etvT300 = sort(rs_etvT300);
boot_etvT400 = sort(rs_etvT400);
boot_etvT500 = sort(rs_etvT500);
boot_etvT700 = sort(rs_etvT700);
boot_etvT850 = sort(rs_etvT850);
boot_etvT925 = sort(rs_etvT925);
boot_etvT1000 = sort(rs_etvT1000);

boot_etvTD150 = sort(rs_etvTD150);
boot_etvTD200 = sort(rs_etvTD200);
boot_etvTD250 = sort(rs_etvTD250);
boot_etvTD300 = sort(rs_etvTD300);
boot_etvTD400 = sort(rs_etvTD400);
boot_etvTD500 = sort(rs_etvTD500);
boot_etvTD700 = sort(rs_etvTD700);
boot_etvTD850 = sort(rs_etvTD850);
boot_etvTD925 = sort(rs_etvTD925);
boot_etvTD1000 = sort(rs_etvTD1000);


boot_etvWIND150 = sort(rs_etvWIND150);
boot_etvWIND200 = sort(rs_etvWIND200);
boot_etvWIND250 = sort(rs_etvWIND250);
boot_etvWIND300 = sort(rs_etvWIND300);
boot_etvWIND400 = sort(rs_etvWIND400);
boot_etvWIND500 = sort(rs_etvWIND500);
boot_etvWIND700 = sort(rs_etvWIND700);
boot_etvWIND850 = sort(rs_etvWIND850);
boot_etvWIND925 = sort(rs_etvWIND925);
boot_etvWIND1000 = sort(rs_etvWIND1000);




siglevel=0.05;  ninetyith = siglevel*resample;
display(['Statistical Significance at bound: ' num2str(1-siglevel)])

sigbound_etvT150 = boot_etvT150(ninetyith)
sigbound_etvT200 = boot_etvT200(ninetyith)
sigbound_etvT250 = boot_etvT250(ninetyith)
sigbound_etvT300 = boot_etvT300(ninetyith)
sigbound_etvT400 = boot_etvT400(ninetyith)
sigbound_etvT500 = boot_etvT500(ninetyith)
sigbound_etvT700 = boot_etvT700(ninetyith)
sigbound_etvT850 = boot_etvT850(ninetyith)
sigbound_etvT925 = boot_etvT925(ninetyith)
sigbound_etvT1000 = boot_etvT1000(ninetyith)


sigbound_etvTD150 = boot_etvTD150(ninetyith)
sigbound_etvTD200 = boot_etvTD200(ninetyith)
sigbound_etvTD250 = boot_etvTD250(ninetyith)
sigbound_etvTD300 = boot_etvTD300(ninetyith)
sigbound_etvTD400 = boot_etvTD400(ninetyith)
sigbound_etvTD500 = boot_etvTD500(ninetyith)
sigbound_etvTD700 = boot_etvTD700(ninetyith)
sigbound_etvTD850 = boot_etvTD850(ninetyith)
sigbound_etvTD925 = boot_etvTD925(ninetyith)
sigbound_etvTD1000 = boot_etvTD1000(ninetyith)


sigbound_etvWIND150 = boot_etvWIND150(ninetyith)
sigbound_etvWIND200 = boot_etvWIND200(ninetyith)
sigbound_etvWIND250 = boot_etvWIND250(ninetyith)
sigbound_etvWIND300 = boot_etvWIND300(ninetyith)
sigbound_etvWIND400 = boot_etvWIND400(ninetyith)
sigbound_etvWIND500 = boot_etvWIND500(ninetyith)
sigbound_etvWIND700 = boot_etvWIND700(ninetyith)
sigbound_etvWIND850 = boot_etvWIND850(ninetyith)
sigbound_etvWIND925 = boot_etvWIND925(ninetyith)
sigbound_etvWIND1000 = boot_etvWIND1000(ninetyith)

%     

sig_etvTheight=[sigbound_etvT1000; sigbound_etvT925; sigbound_etvT850; sigbound_etvT700; sigbound_etvT500;
                sigbound_etvT400; sigbound_etvT300; sigbound_etvT250; sigbound_etvT200; sigbound_etvT150;];
sig_etvTDheight=[sigbound_etvTD1000;sigbound_etvTD925; sigbound_etvTD850; sigbound_etvTD700; sigbound_etvTD500;
                sigbound_etvTD400; sigbound_etvTD300; sigbound_etvTD250; sigbound_etvTD200; sigbound_etvTD150;];
sig_etvWINDheight=[sigbound_etvWIND1000; sigbound_etvWIND925; sigbound_etvWIND850; sigbound_etvWIND700; sigbound_etvWIND500;
                sigbound_etvWIND400; sigbound_etvWIND300; sigbound_etvWIND250; sigbound_etvWIND200; sigbound_etvWIND150;];           
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 



set(groot,'defaultLineLineWidth',2)

height=[1000,925,850,700,500,400,300,250,200,150];
ideal=[1,1,1,1,1,1,1,1,1,1];
zero=[0,0,0,0,0,0,0,0,0,0];

%%% All variables in one plot
% figure
% % subplot(1,2,1)
% h1=plot(D1fixD2fix_etvheight_WIND,height,'b'); hold on; plot(D1varD2fix_etvheight_WIND,height,'b--'); hold on; plot(D1varD2var_etvheight_WIND,height,'b:'); hold on
% h4=plot(sig_etvWINDheight,height,'LineWidth',1,'color','b'); %,[0.4,0.4,0.4])
% hold on
% set(gca,'ydir','reverse')
% h2=plot(D1fixD2fix_etvheight_T,height,'r'); hold on; plot(D1varD2fix_etvheight_T,height,'r--'); hold on; plot(D1varD2var_etvheight_T,height,'r:'); hold on; plot(sig_etvTheight,height,'LineWidth',1,'color','r');
% hold on
% h3=plot(D1fixD2fix_etvheight_TD,height,'g'); hold on; plot(D1varD2fix_etvheight_TD,height,'g--'); hold on; plot(D1varD2var_etvheight_TD,height,'g:'); hold on; plot(sig_etvTDheight,height,'LineWidth',1,'color','g');
% % plot(etvheight_D3WINDFIX,height,'r'); hold on; plot(etvheight_D3WINDVAR1,height,'r--'); hold on; plot(etvheight_D3WINDVAR2,height,'r:') 
% hold on; plot(ideal,height,'k')
% legend([h1 h2 h3 h4], {'Wind','T','TD', '95% sig'},'Location','NorthWest')
% title('Domain 2 Fixed (D1 Fixed:solid; D1 Varied:dashed) vs Varied (D1 Varied:dotted) physics'); hold on; axis tight; xlim([0 3])
% ylabel('Pressure (hPa)')
% xlabel('Error-to-Variance')


% %%% All variables in seperate subplots
% figure
% title('Domain 2 Fixed (D1 Fixed:solid; D1 Varied:dashed) vs Varied (D1 Varied:dotted) physics');
% subplot(1,3,1)
% h1=plot(D1fixD2fix_etvheight_WIND,height,'b'); hold on; plot(D1varD2fix_etvheight_WIND,height,'b--'); hold on; plot(D1varD2var_etvheight_WIND,height,'b:'); hold on
% h4=plot(sig_etvWINDheight,height,'LineWidth',1,'color','b'); %,[0.4,0.4,0.4])
% hold on; axis tight
% hold on; plot(ideal,height,'k')
% set(gca,'ydir','reverse')
% ylabel('Pressure (hPa)')
% xlabel('Error-to-Variance')
% 
% subplot(1,3,2)
% h2=plot(D1fixD2fix_etvheight_T,height,'r'); hold on; plot(D1varD2fix_etvheight_T,height,'r--'); hold on; plot(D1varD2var_etvheight_T,height,'r:'); hold on; plot(sig_etvTheight,height,'LineWidth',1,'color','r');
% hold on
% set(gca,'ydir','reverse'); axis tight
% hold on; plot(ideal,height,'k')
% ylabel('Pressure (hPa)')
% xlabel('Error-to-Variance')
% 
% subplot(1,3,3)
% h3=plot(D1fixD2fix_etvheight_TD,height,'g'); hold on; plot(D1varD2fix_etvheight_TD,height,'g--'); hold on; plot(D1varD2var_etvheight_TD,height,'g:'); hold on; plot(sig_etvTDheight,height,'LineWidth',1,'color','g');
% % plot(etvheight_D3WINDFIX,height,'r'); hold on; plot(etvheight_D3WINDVAR1,height,'r--'); hold on; plot(etvheight_D3WINDVAR2,height,'r:') 
% hold on; plot(ideal,height,'k')
% set(gca,'ydir','reverse')
% % legend([h1 h2 h3 h4], {'Wind','T','TD', '95% sig'},'Location','SouthWest')
% %  hold on; axis tight; xlim([0.5 3])
% ylabel('Pressure (hPa)')
% xlabel('Error-to-Variance')



% 
% height=[1000,925,850,700,500,400,300,250,200,150];
% ideal=[1,1,1,1,1,1,1,1,1,1];
% figure
% % subplot(1,2,1)
% h1=plot(D1fixD2fix_rmseheight_WIND,height,'b'); hold on; plot(D1varD2fix_rmseheight_WIND,height,'b--'); hold on; plot(D1varD2var_rmseheight_WIND,height,'b:')
% hold on
% set(gca,'ydir','reverse')
% h2=plot(D1fixD2fix_rmseheight_T,height,'r'); hold on; plot(D1varD2fix_rmseheight_T,height,'r--'); hold on; plot(D1varD2var_rmseheight_T,height,'r:')
% hold on
% h3=plot(D1fixD2fix_rmseheight_TD,height,'g'); hold on; plot(D1varD2fix_rmseheight_TD,height,'g--'); hold on; plot(D1varD2var_rmseheight_TD,height,'g:')
% legend([h1 h2 h3], {'Wind','T','TD'},'Location','NorthEast')
% title('Domain 2 Fixed (D1 Fixed:solid; D1 Varied:dashed) vs Varied (D1 Varied:dotted) physics'); hold on; axis tight; xlim([0 9])
% ylabel('Pressure (hPa)')
% xlabel('Root Mean Square Error')
% 
% 
% 
% figure
% h1=plot(D1fixD2fix_varianceheight_WIND,height,'b'); hold on; plot(D1varD2fix_varianceheight_WIND,height,'b--'); hold on; plot(D1varD2var_varianceheight_WIND,height,'b:')
% hold on
% set(gca,'ydir','reverse')
% h2=plot(D1fixD2fix_varianceheight_T,height,'r'); hold on; plot(D1varD2fix_varianceheight_T,height,'r--'); hold on; plot(D1varD2var_varianceheight_T,height,'r:')
% hold on
% % h3=plot(D1fixD2fix_varianceheight_TD,height,'g'); hold on; plot(D1varD2fix_varianceheight_TD,height,'g--'); hold on; plot(D1varD2var_varianceheight_TD,height,'g:')
% legend([h1 h2], {'Wind','T'},'Location','North')
% title('Domain 2 Fixed (D1 Fixed:solid; D1 Varied:dashed) vs Varied (D1 Varied:dotted) physics'); hold on; axis tight; xlim([0 18])
% ylabel('Pressure (hPa)')
% xlabel('Variance')


figure('paperunits','in','paperposition',[.5 .5 16 8])

%%%%%%%%%% ALL TIMES
subplot(1,3,1)
height=[1000,925,850,700,500,400,300,250,200,150];
zero=[0,0,0,0,0,0,0,0,0,0];
ideal=[1,1,1,1,1,1,1,1,1,1];
h1=plot(D1fixD2fix_etvheight_WIND,height,'b'); hold on; plot(D1varD2fix_etvheight_WIND,height,'b--');set(gca,'FontSize',16); hold on; %plot(D1varD2var_etvheight_WIND,height,'b:'); hold on
h4=plot(sig_etvWINDheight,height,'LineWidth',1,'color','b'); %,[0.4,0.4,0.4])
hold on
set(gca,'ydir','reverse')
h2=plot(D1fixD2fix_etvheight_T,height,'r'); hold on; plot(D1varD2fix_etvheight_T,height,'r--'); hold on; plot(sig_etvTheight,height,'LineWidth',1,'color','r'); %hold on; plot(D1varD2var_etvheight_T,height,'r:');
hold on
h3=plot(D1fixD2fix_etvheight_TD,height,'g'); hold on; plot(D1varD2fix_etvheight_TD,height,'g--'); hold on; plot(sig_etvTDheight,height,'LineWidth',1,'color','g'); %hold on; plot(D1varD2var_etvheight_TD,height,'g:');
% plot(etvheight_D3WINDFIX,height,'r'); hold on; plot(etvheight_D3WINDVAR1,height,'r--'); hold on; plot(etvheight_D3WINDVAR2,height,'r:') 
hold on; plot(ideal,height,'k')
% legend([h1 h2 h3 h4], {'Wind','T','TD', '95% sig'},'Location','East')
title('All Forecast hours 0-48'); hold on; axis tight; xlim([0 4])
ylabel('Pressure (hPa)','fontsize', 16)
xlabel('Error-to-Variance','fontsize', 16)


%%%%%%%%%% 12-24 hours
subplot(1,3,2)
height=[1000,925,850,700,500,400,300,250,200,150];
ideal=[1,1,1,1,1,1,1,1,1,1];
h1=plot(D1fixD2fix_etvheight_WIND,height,'b'); hold on; plot(D1varD2fix_etvheight_WIND,height,'b--');set(gca,'FontSize',16); hold on; %plot(D1varD2var_etvheight_WIND,height,'b:'); hold on
h4=plot(sig_etvWINDheight,height,'LineWidth',1,'color','b'); %,[0.4,0.4,0.4])
hold on
set(gca,'ydir','reverse')
h2=plot(D1fixD2fix_etvheight_T,height,'r'); hold on; plot(D1varD2fix_etvheight_T,height,'r--'); hold on; plot(sig_etvTheight,height,'LineWidth',1,'color','r');hold on; %plot(D1varD2var_etvheight_T,height,'r:');
hold on
h3=plot(D1fixD2fix_etvheight_TD,height,'g'); hold on; plot(D1varD2fix_etvheight_TD,height,'g--'); hold on; plot(sig_etvTDheight,height,'LineWidth',1,'color','g'); hold on; %plot(D1varD2var_etvheight_TD,height,'g:');
% plot(etvheight_D3WINDFIX,height,'r'); hold on; plot(etvheight_D3WINDVAR1,height,'r--'); hold on; plot(etvheight_D3WINDVAR2,height,'r:') 
hold on; plot(ideal,height,'k')
% legend([h1 h2 h3 h4], {'Wind','T','TD', '95% sig'},'Location','NorthWest')
title('Forecast hours 12-24'); hold on; axis tight; xlim([0 4])
ylabel('Pressure (hPa)','fontsize', 16)
xlabel('Error-to-Variance','fontsize', 16)



%%%%%%%%%% 36-48 hours
subplot(1,3,3)
height=[1000,925,850,700,500,400,300,250,200,150];
ideal=[1,1,1,1,1,1,1,1,1,1];
h1=plot(D1fixD2fix_etvheight_WIND,height,'b'); hold on; plot(D1varD2fix_etvheight_WIND,height,'b--');set(gca,'FontSize',16); hold on; %plot(D1varD2var_etvheight_WIND,height,'b:'); hold on
h4=plot(sig_etvWINDheight,height,'LineWidth',1,'color','b'); %,[0.4,0.4,0.4])
hold on
set(gca,'ydir','reverse')
h2=plot(D1fixD2fix_etvheight_T,height,'r'); hold on; plot(D1varD2fix_etvheight_T,height,'r--'); hold on; plot(sig_etvTheight,height,'LineWidth',1,'color','r'); hold on; %plot(D1varD2var_etvheight_T,height,'r:');
hold on
h3=plot(D1fixD2fix_etvheight_TD,height,'g'); hold on; plot(D1varD2fix_etvheight_TD,height,'g--'); hold on; plot(sig_etvTDheight,height,'LineWidth',1,'color','g'); %hold on; plot(D1varD2var_etvheight_TD,height,'g:');
% plot(etvheight_D3WINDFIX,height,'r'); hold on; plot(etvheight_D3WINDVAR1,height,'r--'); hold on; plot(etvheight_D3WINDVAR2,height,'r:') 
hold on; plot(ideal,height,'k')
legend([h1 h2 h3 h4], {'Wind','T','TD', '95% sig'},'Location','East', 'fontsize', 15)
title('Forecast hours 36-48'); hold on; axis tight; hold on; xlim([0 4])
ylabel('Pressure (hPa)','fontsize', 16)
xlabel('Error-to-Variance','fontsize', 16)



set(gcf,'units','normalized','outerposition',[0 0 1 1])
hold on
annotation('textbox', [0 0.9 1 0.1],'FontSize',20,'String', 'Domain 2 Fixed (D1 Fixed:solid; D1 Varied:dashed) physics','EdgeColor', 'none','HorizontalAlignment', 'center')

% print('D2_fix_var_etv_disagg_raob_times.png','-dpng','-r400')
figure('paperunits','in','paperposition',[.5 .5 16 8])




%%%%%%%%%% ALL TIMES
figure
height=[1000,925,850,700,500,400,300,250,200,150];
ideal=[1,1,1,1,1,1,1,1,1,1];
h1=plot(D1fixD2fix_etvheight_WIND,height,'b'); hold on; plot(D1varD2fix_etvheight_WIND,height,'b--');set(gca,'FontSize',16); hold on; plot(D1varD2var_etvheight_WIND,height,'b:'); hold on
h4=plot(sig_etvWINDheight,height,'LineWidth',1,'color','b'); %,[0.4,0.4,0.4])
hold on
set(gca,'ydir','reverse')
h2=plot(D1fixD2fix_etvheight_T,height,'r'); hold on; plot(D1varD2fix_etvheight_T,height,'r--'); hold on; plot(sig_etvTheight,height,'LineWidth',1,'color','r'); hold on; plot(D1varD2var_etvheight_T,height,'r:');
hold on
h3=plot(D1fixD2fix_etvheight_TD,height,'g'); hold on; plot(D1varD2fix_etvheight_TD,height,'g--'); hold on; plot(sig_etvTDheight,height,'LineWidth',1,'color','g'); hold on; plot(D1varD2var_etvheight_TD,height,'g:');
% plot(etvheight_D3WINDFIX,height,'r'); hold on; plot(etvheight_D3WINDVAR1,height,'r--'); hold on; plot(etvheight_D3WINDVAR2,height,'r:') 
hold on; plot(ideal,height,'k')
% legend([h1 h2 h3 h4], {'Wind','T','TD', '95% sig'},'Location','East')
axis tight; xlim([0 3.5])
ylabel('Pressure (hPa)','fontsize', 16)
xlabel('Error-to-Variance','fontsize', 16)
set(gcf,'units','normalized','outerposition',[0 0 1 1])






%%%%%% including PBL and MP only varied. plot RMSE, VARIANCE, AND BIAS
figure('paperunits','in','paperposition',[.5 .5 16 8])
height=[1000,925,850,700,500,400,300,250,200,150];
ideal=[1,1,1,1,1,1,1,1,1,1];
h1=plot(D1varD2var_etvheight_WIND,height,'b'); hold on; plot(D2PBL_etvheight_WIND,height,'b--'); hold on; plot(D2MP_etvheight_WIND,height,'b:'); hold on
set(gca,'FontSize',16); hold on
% h4=plot(sig_etvWINDheight,height,'LineWidth',1,'color','b'); %,[0.4,0.4,0.4])
hold on
set(gca,'ydir','reverse')
h2=plot(D1varD2var_etvheight_T,height,'r'); hold on; plot(D2PBL_etvheight_T,height,'r--'); hold on; plot(D2MP_etvheight_T,height,'r:');
hold on
h3=plot(D1varD2var_etvheight_TD,height,'g'); hold on; plot(D2PBL_etvheight_TD,height,'g--'); hold on; plot(D2MP_etvheight_TD,height,'g:');
% plot(etvheight_D3WINDFIX,height,'r'); hold on; plot(etvheight_D3WINDVAR1,height,'r--'); hold on; plot(etvheight_D3WINDVAR2,height,'r:') 
hold on; plot(ideal,height,'k')
% legend([h1 h2 h3 h4], {'Wind','T','TD', '95% sig'},'Location','East')
axis tight; xlim([0 3.5])
ylabel('Pressure (hPa)','fontsize', 16)
xlabel('Error-to-Variance','fontsize', 16)
set(gcf,'units','normalized','outerposition',[0 0 1 1])






set(groot,'defaultLineLineWidth',2)
height=[1000,925,850,700,500,400,300,250,200,150];
zero=[0,0,0,0,0,0,0,0,0,0];
figure('paperunits','in','paperposition',[.5 .5 16 8])
subplot(1,3,1)
h1=plot(D1varD2var_rmseheight_WIND,height,'b'); hold on; plot(D2PBL_rmseheight_WIND,height,'b--'); hold on; plot(D2MP_rmseheight_WIND,height,'b:');
hold on
set(gca,'FontSize',16);
set(gca,'ydir','reverse')
h2=plot(D1varD2var_rmseheight_T,height,'r'); hold on; plot(D2PBL_rmseheight_T,height,'r--'); hold on;plot(D2MP_rmseheight_T,height,'r:'); 
hold on
% h3=plot(D1varD2var_rmseheight_TD,height,'g'); hold on; plot(D2PBL_rmseheight_TD,height,'g--'); hold on; plot(D2MP_rmseheight_TD,height,'g:');
% legend([h1 h2 h3], {'Wind','T', 'TD'},'Location','NorthEast')
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics'); hold on; 
axis tight; %xlim([0 5])
ylabel('Pressure (hPa)', 'fontsize', 15)
xlabel('RMSE', 'fontsize', 15)
axis tight

% figure
subplot(1,3,2)
h1=plot(D1varD2var_varianceheight_WIND,height,'b'); hold on; plot(D2PBL_varianceheight_WIND,height,'b--'); hold on; plot(D2MP_varianceheight_WIND,height,'b:');
hold on
set(gca,'FontSize',16);
set(gca,'ydir','reverse')
h2=plot(D1varD2var_varianceheight_T,height,'r'); hold on; plot(D2PBL_varianceheight_T,height,'r--'); hold on;plot(D2MP_varianceheight_T,height,'r:'); 
hold on
% h3=plot(D1varD2var_varianceheight_TD,height,'g'); hold on; plot(D2PBL_varianceheight_TD,height,'g--'); hold on; plot(D2MP_varianceheight_TD,height,'g:');
% legend([h1 h2 h3], {'Wind','T', 'TD'},'Location','NorthEast')
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics'); hold on; 
axis tight; %xlim([0 20])
ylabel('Pressure (hPa)', 'fontsize', 15)
xlabel('Variance', 'fontsize', 15)
axis tight

% figure
subplot(1,3,3)
h1=plot(D1varD2var_biasheight_WIND,height,'b'); hold on; plot(D2PBL_biasheight_WIND,height,'b--'); hold on; plot(D2MP_biasheight_WIND,height,'b:');
hold on
set(gca,'FontSize',16);
set(gca,'ydir','reverse')
h2=plot(D1varD2var_biasheight_T,height,'r'); hold on; plot(D2PBL_biasheight_T,height,'r--'); hold on;plot(D2MP_biasheight_T,height,'r:'); 
hold on
% h3=plot(D1varD2var_biasheight_TD,height,'g'); hold on; plot(D2PBL_biasheight_TD,height,'g--'); hold on; plot(D2MP_biasheight_TD,height,'g:');
hold on; plot(zero,height,'k')
legend([h1 h2], {'Wind','T'},'Location','SouthEast', 'fontsize', 15) %h3 , 'TD'
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics'); hold on; axis tight; %xlim([0 7])
ylabel('Pressure (hPa)', 'fontsize', 15)
xlabel('Bias', 'fontsize', 15)
axis tight

set(gcf,'units','normalized','outerposition',[0 0 1 1])











%%%%%% plot RMSE, VARIANCE, AND BIAS for fixfix varfix varvar
set(groot,'defaultLineLineWidth',2)
height=[1000,925,850,700,500,400,300,250,200,150];
zero=[0,0,0,0,0,0,0,0,0,0];
figure('paperunits','in','paperposition',[.5 .5 16 8])
subplot(1,3,1)
h1=plot(D1fixD2fix_rmseheight_WIND,height,'b'); hold on; plot(D1varD2fix_rmseheight_WIND,height,'b--'); hold on; plot(D1varD2var_rmseheight_WIND,height,'b:');
hold on
set(gca,'FontSize',16);
set(gca,'ydir','reverse')
h2=plot(D1fixD2fix_rmseheight_T,height,'r'); hold on; plot(D1varD2fix_rmseheight_T,height,'r--'); hold on;plot(D1varD2var_rmseheight_T,height,'r:'); 
hold on
h3=plot(D1fixD2fix_rmseheight_TD,height,'g'); hold on; plot(D1varD2fix_rmseheight_TD,height,'g--'); hold on; plot(D1varD2var_rmseheight_TD,height,'g:');
% legend([h1 h2 h3], {'Wind','T', 'TD'},'Location','NorthEast')
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics'); hold on; 
axis tight; %xlim([0 5])
ylabel('Pressure (hPa)', 'fontsize', 15)
xlabel('RMSE', 'fontsize', 15)
axis tight

% figure
subplot(1,3,2)
h1=plot(D1fixD2fix_varianceheight_WIND,height,'b'); hold on; plot(D1varD2fix_varianceheight_WIND,height,'b--'); hold on; plot(D1varD2var_varianceheight_WIND,height,'b:');
hold on
set(gca,'FontSize',16);
set(gca,'ydir','reverse')
h2=plot(D1fixD2fix_varianceheight_T,height,'r'); hold on; plot(D1varD2fix_varianceheight_T,height,'r--'); hold on;plot(D1varD2var_varianceheight_T,height,'r:'); 
hold on
h3=plot(D1fixD2fix_varianceheight_TD,height,'g'); hold on; plot(D1varD2fix_varianceheight_TD,height,'g--'); hold on; plot(D1varD2var_varianceheight_TD,height,'g:');
% legend([h1 h2 h3], {'Wind','T', 'TD'},'Location','NorthEast')
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics'); hold on; 
axis tight; %xlim([0 20])
ylabel('Pressure (hPa)', 'fontsize', 15)
xlabel('Variance', 'fontsize', 15)
axis tight

% figure
subplot(1,3,3)
h1=plot(D1fixD2fix_biasheight_WIND,height,'b'); hold on; plot(D1varD2fix_biasheight_WIND,height,'b--'); hold on; plot(D1varD2var_biasheight_WIND,height,'b:');
hold on
set(gca,'FontSize',16);
set(gca,'ydir','reverse')
h2=plot(D1fixD2fix_biasheight_T,height,'r'); hold on; plot(D1varD2fix_biasheight_T,height,'r--'); hold on;plot(D1varD2var_biasheight_T,height,'r:'); 
hold on
h3=plot(D1fixD2fix_biasheight_TD,height,'g'); hold on; plot(D1varD2fix_biasheight_TD,height,'g--'); hold on; plot(D1varD2var_biasheight_TD,height,'g:');
hold on; plot(zero,height,'k')
legend([h1 h2 h3], {'Wind','T', 'TD'},'Location','SouthEast', 'fontsize', 15)
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics'); hold on; axis tight; %xlim([0 7])
ylabel('Pressure (hPa)', 'fontsize', 15)
xlabel('Bias', 'fontsize', 15)
axis tight

set(gcf,'units','normalized','outerposition',[0 0 1 1])














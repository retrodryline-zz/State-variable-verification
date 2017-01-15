clearvars -except D1* D2* sigbound*
% clearvars -except *time_mem*
 
spatial_agg = 0; 
outofpoly=0 %1 does outside of poly, 0 inside


%domain 2 west edge
% poly_y = [27.75,41.9,41.9,27.75]; %sw,nw,ne,se corners 
% poly_x = [-104.25,-105.6,-102.5,-102.55]; %sw,nw,ne,se corners



%domain 1 LBCs (400 km)
% poly_y = [26.05, 48.88, 51.9, 27.07]; %sw,nw,ne,se corners 
% poly_x = [-117.15,-127.18,-84.45,-85.4]; %sw,nw,ne,se corners

%domain 1 LBCs (600 km)
% poly_y = [27.02, 46.0, 48.05, 28.47];% poly_y = [27.6, 47.55, 49.8, 28.55]; %sw,nw,ne,se corners 
% poly_x = [-115.97,-122.69,-85.4,-89.4];% poly_x = [-115.85,-123.55,-85.95,-88.2]; %sw,nw,ne,se corners

%domain 1 LBCs (800 km)
% poly_y = [29.95, 47.0, 47.16, 30.4]; %sw,nw,ne,se corners 
% poly_x = [-115.95,-119.68,-89.45,-89.75]; %sw,nw,ne,se corners

%domain 1 LBCs (600 km) JUST WESTERN SIDE (assuming in poly)
% poly_y = [17.9, 50.0, 55.0, 17.8];% poly_y = [19.8, 54.45, 57.0, 18.0, 19.8]; %sw,nw,ne,se corners 
% poly_x = [-122.0, -139.0,-125.5,-112.4];% poly_x = [-113.65,-126.3,-69.0,-73.0, -113.65]; %sw,nw,ne,se corners


%domain 1 LBCs (600 km) JUST EASTERN (assuming in  poly)
% poly_y = [21.3, 54.0, 52.7, 20.5];% poly_y = [21.65, 56.5, 54.9, 20.5, 21.65]; %sw,nw,ne,se corners 
% poly_x = [-90.9, -84.3, -71.0, -79.8];% poly_x = [-89.0, -85.6, -70.7, -79.8, -89.0]; %sw,nw,ne,se corners


%domain 2 LBCS (300 km)
% poly_y = [31.3, 38.8, 38.8, 31.18]; %sw,nw,ne,se corners 
% poly_x = [-100.7, -101.16, -95.13, -95.52]; %sw,nw,ne,se corners

%domain 2 LBCS (200 km)
% poly_y = [30.24, 39.67, 39.7, 30.27]; %sw,nw,ne,se corners 
% poly_x = [-101.73, -102.37, -93.9, -94.56]; %sw,nw,ne,se corners

%domain 2 LBCS (100 km)
% poly_y = [29.34, 40.6, 40.56, 29.4]; %sw,nw,ne,se corners 
% poly_x = [-102.68, -103.6, -92.7, -93.58]; %sw,nw,ne,se corners


start=pwd;
 
num_mems=40; %for indexing
num_times=49; %forecast length + 1

sigma_T2=1.5; %default 2.5
sigma_TD2=1.5; %default 2.5 %RH=0.2 error
sigma_WIND10=1.5; %default 3.5

% D1var_errtovarT2_time=errtovarT2_time;
% D1var_errtovarTD2_time=errtovarTD2_time;
% D1var_errtovarWIND10_time=errtovarWIND10_time;
% D1var_varT2_time=varT2_time;
% D1var_varTD2_time=varTD2_time;
% D1var_varWIND10_time=varWIND10_time;
% D1var_RMST2_time=RMST2_time;
% D1var_RMSTD2_time=RMSTD2_time;
% D1var_RMSWIND10_time=RMSWIND10_time;
% D1var_BIAST2_time=biasT2_time;
% D1var_BIASTD2_time=biasTD2_time;
% D1var_BIASWIND10_time=biasWIND10_time;

% D1PBL_errtovarT2_time=errtovarT2_time;
% D1PBL_errtovarTD2_time=errtovarTD2_time;
% D1PBL_errtovarWIND10_time=errtovarWIND10_time;
% D1PBL_varT2_time=varT2_time;
% D1PBL_varTD2_time=varTD2_time;
% D1PBL_varWIND10_time=varWIND10_time;
% D1PBL_RMST2_time=RMST2_time;
% D1PBL_RMSTD2_time=RMSTD2_time;
% D1PBL_RMSWIND10_time=RMSWIND10_time;
% D1PBL_BIAST2_time=biasT2_time;
% D1PBL_BIASTD2_time=biasTD2_time;
% D1PBL_BIASWIND10_time=biasWIND10_time;

% 
% 
% D2MP_errtovarT2_time=errtovarT2_time;
% D2MP_errtovarTD2_time=errtovarTD2_time;
% D2MP_errtovarWIND10_time=errtovarWIND10_time;
% D2MP_varT2_time=varT2_time;
% D2MP_varTD2_time=varTD2_time;
% D2MP_varWIND10_time=varWIND10_time;
% D2MP_RMST2_time=RMST2_time;
% D2MP_RMSTD2_time=RMSTD2_time;
% D2MP_RMSWIND10_time=RMSWIND10_time;
% D2MP_BIAST2_time=biasT2_time;
% D2MP_BIASTD2_time=biasTD2_time;
% D2MP_BIASWIND10_time=biasWIND10_time;
% D2MP_RMST2_time_mem=RMST2_mem;
% D2MP_RMSTD2_time_mem=RMSTD2_mem;
% D2MP_RMSWIND10_time_mem=RMSWIND10_mem;
% D2MP_ERRT2_time_mem=ERRT2_mem;
% D2MP_ERRTD2_time_mem=ERRTD2_mem;
% D2MP_ERRWIND10_time_mem=ERRWIND10_mem;
% D2MP_MEANDIFFT2_time_mem=MEANDIFFT2_mem;
% D2MP_MEANDIFFTD2_time_mem=MEANDIFFTD2_mem;
% D2MP_MEANDIFFWIND10_time_mem=MEANDIFFWIND10_mem;





% D1varD2var_errtovarT2_time=errtovarT2_time;
% D1varD2var_errtovarTD2_time=errtovarTD2_time;
% D1varD2var_errtovarWIND10_time=errtovarWIND10_time;
% D1varD2var_varT2_time=varT2_time;
% D1varD2var_varTD2_time=varTD2_time;
% D1varD2var_varWIND10_time=varWIND10_time;
% D1varD2var_RMST2_time=RMST2_time;
% D1varD2var_RMSTD2_time=RMSTD2_time;
% D1varD2var_RMSWIND10_time=RMSWIND10_time;
% D1varD2var_BIAST2_time=biasT2_time;
% D1varD2var_BIASTD2_time=biasTD2_time;
% D1varD2var_BIASWIND10_time=biasWIND10_time;
% D1varD2var_RMST2_time_mem=RMST2_mem;
% D1varD2var_RMSTD2_time_mem=RMSTD2_mem;
% D1varD2var_RMSWIND10_time_mem=RMSWIND10_mem;

% D1fixD2fix_ERRT2_time_mem=ERRT2_mem;
% D1fixD2fix_ERRTD2_time_mem=ERRTD2_mem;
% D1fixD2fix_ERRWIND10_time_mem=ERRWIND10_mem;

% D1varD2var_MEANDIFFT2_time_mem=MEANDIFFT2_mem;
% D1varD2var_MEANDIFFTD2_time_mem=MEANDIFFTD2_mem;
% D1varD2var_MEANDIFFWIND10_time_mem=MEANDIFFWIND10_mem;


% % % % D1fixD2fix_etvT2 = errtovarT2_case;
% % % % D1fixD2fix_etvTD2 = errtovarTD2_case;
% % % % D1fixD2fix_etvWIND10 = errtovarWIND10_case;

type='var'; %set this to var for d1var to avoid resampling sig stats
d1fov='var';
d1fov_dir='MPvar';

domain_list=[
%         'D1'
        'D2'
        ]; %don't do both domains

list=[
    '2015050600'
    '2015050612'
    '2015050700'
    '2015050712'
    '2015050800'
    '2015050812'
    '2015050900'
    '2015050912'   
    ];    %050700 had bad iowa TD2 ob
%%%%%%%%%%%%%%%%%%%% ALLOCATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RMST2_timemem=zeros(num_mems,num_times,size(list,1)); RMSTD2_timemem=zeros(num_mems,num_times,size(list,1)); RMSWIND10_timemem=zeros(num_mems,num_times,size(list,1));
RMST2_mem=zeros(num_mems,num_times); RMSTD2_mem=zeros(num_mems,num_times); RMSWIND10_mem=zeros(num_mems,num_times);

ERRT2_timemem=zeros(num_mems,num_times,size(list,1)); ERRTD2_timemem=zeros(num_mems,num_times,size(list,1)); ERRWIND10_timemem=zeros(num_mems,num_times,size(list,1));
ERRT2_mem=zeros(num_mems,num_times); ERRTD2_mem=zeros(num_mems,num_times); ERRWIND10_mem=zeros(num_mems,num_times);

MEANDIFFT2_timemem=zeros(num_mems,num_times,size(list,1)); MEANDIFFTD2_timemem=zeros(num_mems,num_times,size(list,1)); MEANDIFFWIND10_timemem=zeros(num_mems,num_times,size(list,1));
MEANDIFFT2_mem=zeros(num_mems,num_times); MEANDIFFTD2_mem=zeros(num_mems,num_times); MEANDIFFWIND10_mem=zeros(num_mems,num_times);



biasT2=zeros(num_times,size(list,1)); biasTD2=zeros(num_times,size(list,1)); biasWIND10=zeros(num_times,size(list,1));
RMST2=zeros(num_times,size(list,1)); RMSTD2=zeros(num_times,size(list,1)); RMSWIND10=zeros(num_times,size(list,1));
varT2=zeros(num_times,size(list,1)); varTD2=zeros(num_times,size(list,1)); varWIND10=zeros(num_times,size(list,1));
errtovarT2=zeros(num_times,size(list,1)); errtovarTD2=zeros(num_times,size(list,1)); errtovarWIND10=zeros(num_times,size(list,1));


obs_count_T2_rel=zeros(num_times,size(list,1)); obs_count_TD2_rel=zeros(num_times,size(list,1)); obs_count_WIND10_rel=zeros(num_times,size(list,1));

errtovarT2_time=zeros(num_times); varT2_time=zeros(num_times); RMST2_time=zeros(num_times); biasT2_time=zeros(num_times);
errtovarTD2_time=zeros(num_times); varTD2_time=zeros(num_times); RMSTD2_time=zeros(num_times); biasTD2_time=zeros(num_times);
errtovarWIND10_time=zeros(num_times); varWIND10_time=zeros(num_times); RMSWIND10_time=zeros(num_times); biasWIND10_time=zeros(num_times);
errtovarT2_case=zeros(size(list,1),1); errtovarTD2_case=zeros(size(list,1),1); errtovarWIND10_case=zeros(size(list,1),1);

rank_sum_T2=zeros(num_mems+1,num_times); rank_sum_TD2=zeros(num_mems+1,num_times); rank_sum_WIND10=zeros(num_mems+1,num_times);    

obs_count_T2=zeros(num_times,size(list,1));
obs_count_TD2=zeros(num_times,size(list,1));
obs_count_WIND10=zeros(num_times,size(list,1));

if strcmp(type,'fix')==1 && strcmp(d1fov,'fix')==1
sigbound_etvT2=zeros(num_times,1); sigbound_etvTD2=zeros(num_times,1); sigbound_etvWIND10=zeros(num_times,1);
end



for y=1:size(domain_list,1)
    domain=domain_list(y,:);
cd([start '/D1' d1fov_dir])    %d1fov
  display(['Working in directory: ' d1fov_dir]) %d1fov
for t=1:num_times 
    
    %summary of stats for all cases for resampling
    RMST2_all_sumcase = [];
    tot_varT2_sumcase = [];
    RMSTD2_all_sumcase = [];
    tot_varTD2_sumcase = [];
    RMSWIND10_all_sumcase = [];
    tot_varWIND10_sumcase = []; 
    
for z=1:size(list,1);
        datestring=list(z,:);    
    if strcmp(domain,'D1')==1
        %doing domain1
        fname=['verifsfc' domain '_' datestring '_f' num2str(t-1)]
    else
        %doing domain2
        fname=['verifsfc' domain type '_' datestring '_f' num2str(t-1)]
    end
    ncid=netcdf.open(fname,'NOWRITE');
    
    
    
    %%% Get LAT and LON info
    varid=netcdf.inqVarID(ncid,'LAT');
    LAT=netcdf.getVar(ncid,varid,'double');
    
    varid=netcdf.inqVarID(ncid,'LON');
    LON=netcdf.getVar(ncid,varid,'double');

    
    
    
    %%% T2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    varid=netcdf.inqVarID(ncid,'T2');
    T2=netcdf.getVar(ncid,varid,'double');

    T2(T2<255)=nan;
    T2(T2>317)=nan;

    [row,col]=find(isnan(T2));
    blahinds=unique(row);
    T2(blahinds,:)=[];
%     T2(isnan(T2(:,1)), :)=[]; %2162
    latT2=LAT;    latT2(blahinds,:)=[];
    lonT2=LON;    lonT2(blahinds,:)=[];
    
    [dummy uniqueT2]=unique([T2(:,1)]);
    T2=T2(uniqueT2,:);
    latT2=latT2(uniqueT2); lonT2=lonT2(uniqueT2);
    

    modobsT2=zeros(size(T2,1),2); %Allocate
    
    for i=1:size(T2,1)
    modobsT2(i,1)=nanmean(T2(i,1:num_mems)); %take ens mean
    modobsT2(i,2)=T2(i,num_mems+1); %no mean just raw obs value
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% TD2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    varid=netcdf.inqVarID(ncid,'TD2');
    TD2=netcdf.getVar(ncid,varid,'double');

    TD2(TD2<250)=nan;
    TD2(TD2>305)=nan;
    
    [row,col]=find(isnan(TD2));
    blahinds=unique(row);
    TD2(blahinds,:)=[];    
%     TD2(isnan(TD2(:,1)), :)=[]; %2162
    latTD2=LAT;    latTD2(blahinds,:)=[];
    lonTD2=LON;    lonTD2(blahinds,:)=[];
    
    [dummy uniqueTD2]=unique([TD2(:,1)]);
    TD2=TD2(uniqueTD2,:);
    latTD2=latTD2(uniqueTD2); lonTD2=lonTD2(uniqueTD2);
    
    modobsTD2=zeros(size(TD2,1),2); %Allocate
    
    for i=1:size(TD2,1)
    modobsTD2(i,1)=nanmean(TD2(i,1:num_mems)); %take ens mean
    modobsTD2(i,2)=TD2(i,num_mems+1); %no mean just raw obs value
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%% U10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    varid=netcdf.inqVarID(ncid,'U10');
    U10=netcdf.getVar(ncid,varid,'double');
    U10(U10<-40 | U10==0)=nan;
    U10(U10>40)=nan;
    
    varid=netcdf.inqVarID(ncid,'V10');
    V10=netcdf.getVar(ncid,varid,'double');
    netcdf.close(ncid)
    V10(V10<-40 | V10==0)=nan;
    V10(V10>40)=nan;
    
    [rowu,colu]=find(isnan(U10));
    [rowv,colv]=find(isnan(V10));
    rowcomb = rowu; rowcomb(size(rowcomb,1)+1:size(rowcomb,1)+size(rowv,1))=rowv;
    blahinds=unique(rowcomb);
    
%     blahinds=unique(row);
    U10(blahinds,:)=[];
    V10(blahinds,:)=[];
%     U10(isnan(U10(:,1)), :)=[];
    latU10=LAT;    latU10(blahinds,:)=[];
    lonU10=LON;    lonU10(blahinds,:)=[];
    latV10=LAT;    latV10(blahinds,:)=[];
    lonV10=LON;    lonV10(blahinds,:)=[];

    modobsU10=zeros(size(U10,1),2); %Allocate
    modobsV10=zeros(size(V10,1),2); %Allocate
    
    for i=1:size(U10,1)
        modobsU10(i,1)=nanmean(U10(i,1:num_mems)); %take ens mean
        modobsU10(i,2)=U10(i,num_mems+1); %no mean just raw obs value
    end
    
    for i=1:size(V10,1)
        modobsV10(i,1)=nanmean(V10(i,1:num_mems)); %take ens mean
        modobsV10(i,2)=V10(i,num_mems+1); %no mean just raw obs value
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% V10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     blahinds=unique(row);               
%     V10(isnan(V10(:,1)), :)=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear ncid


    %%% Calculate magnitude WS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WIND10=sqrt(U10.^2 + V10.^2);

    % remove any observed 0 wind
%     WIND10(WIND10(:,num_mems+1)==0)=nan;
%     WIND10(isnan(WIND10(:,num_mems+1)),:)=[];
%     latU10(blahinds,:)=[]; lonU10(blahinds,:)=[];
    
    %remove any duplicates (easy to do member1)
    [dummy uniqueWIND10]=unique([WIND10(:,1)]);
    WIND10=WIND10(uniqueWIND10,:);
    latU10=latU10(uniqueWIND10); lonU10=lonU10(uniqueWIND10);
  
    
    modobsWIND10=zeros(size(WIND10,1),2); %Allocate
    for i=1:size(WIND10,1)
        modobsWIND10(i,1)=nanmean(WIND10(i,1:num_mems)); %take ens mean
        modobsWIND10(i,2)=WIND10(i,num_mems+1); %no mean just raw obs value
    end
%     modobsWIND10=sqrt(modobsU10.^2 + modobsV10.^2);
%     modobsWIND10(modobsWIND10==0)=nan; %neglect 0 WS obs
%     modobsWIND10(isnan(modobsWIND10(:,2)), :)=[];
    
    latWIND10 = latU10; lonWIND10 = lonU10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    
    %%%%%%%%% SPATIAL DISAGREGATION %%%%%%%%%%%%%%
    if spatial_agg == 1
        display('Doing spatial aggregation based upon input polygon points')
        
        if outofpoly==1 %%% only obtain verification OUTSIDE box (near LBC)
           
            [in,on] = inpolygon(lonT2,latT2,poly_x,poly_y);
            T2 = T2(~in,:);  modobsT2 = modobsT2(~in,:);  latT2 = latT2(~in);  lonT2 = lonT2(~in);
            clearvars in on

%             latlim=[22.0 53.0];
%             lonlim=[-133.0 -77.0];
%             hnd=figure;
%             ax=worldmap(latlim,lonlim);
%             latlim=getm(ax,'MapLatLimit');
%             lonlim=getm(ax,'MapLonLimit');
%             coast=load('coast');
%             geoshow(ax, latlim, lonlim,'DisplayType', 'polygon', 'FaceColor', 'none')
%             states=shaperead('usastatehi','UseGeocoords',true,'BoundingBox',[lonlim',latlim']);
%             geoshow(ax,states,'FaceColor','none')
%             close(hnd);
%             
%             
%             scrsz = get(0,'ScreenSize');
%             figure('Position',[1 scrsz(4)/1.5 scrsz(3)/1.5 scrsz(4)/1.5])
%             hold on
%             xlim(lonlim)
%             ylim(latlim)
%             plot(lonT2(in),latT2(in),'r+')
%             plot(lonT2(~in),latT2(~in),'bo')
%             geoshow(coast,'color','k')
%             for is=1:size(states,1)
%                 plot(states(is).Lon,states(is).Lat,'k');    
%             end

            [in,on] = inpolygon(lonTD2,latTD2,poly_x,poly_y);
            TD2 = TD2(~in,:);    modobsTD2 = modobsTD2(~in,:);  latTD2 = latTD2(~in);  lonTD2 = lonTD2(~in);
            clearvars in on

            [in,on] = inpolygon(lonWIND10,latWIND10,poly_x,poly_y);
            WIND10 = WIND10(~in,:);  modobsWIND10 = modobsWIND10(~in,:);  latWIND10 = latWIND10(~in);  lonWIND10 = lonWIND10(~in);
            clearvars in on
            
        else %%%only obtain verification WITHIN box
            
            [in,on] = inpolygon(lonT2,latT2,poly_x,poly_y);
            T2 = T2(in,:);  modobsT2 = modobsT2(in,:);  latT2 = latT2(in);  lonT2 = lonT2(in);
            clearvars in on
            
            [in,on] = inpolygon(lonTD2,latTD2,poly_x,poly_y);
            TD2 = TD2(in,:);    modobsTD2 = modobsTD2(in,:);  latTD2 = latTD2(in);  lonTD2 = lonTD2(in);
            clearvars in on

            [in,on] = inpolygon(lonWIND10,latWIND10,poly_x,poly_y);
            WIND10 = WIND10(in,:);  modobsWIND10 = modobsWIND10(in,:);  latWIND10 = latWIND10(in);  lonWIND10 = lonWIND10(in);
            clearvars in on
            
        end
    else
        display('Doing verification on entire domain')
    end
    
    obs_count_T2(t,z) = size(T2,1);
    obs_count_TD2(t,z) = size(TD2,1);
    obs_count_WIND10(t,z) = size(WIND10,1);

    
    
    
    
%%%%%%%%%%%%%%%% STATISTICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% RMSE individual mems %%%
    for e=1:num_mems       
        RMST2_timemem(e,t,z)=mean(sqrt((T2(:,e)-T2(:,num_mems+1)).^2));
        RMSTD2_timemem(e,t,z)=mean(sqrt((TD2(:,e)-TD2(:,num_mems+1)).^2));  
        RMSWIND10_timemem(e,t,z)=mean(sqrt((WIND10(:,e)-WIND10(:,num_mems+1)).^2));

        ERRT2_timemem(e,t,z)=mean(T2(:,e)-T2(:,num_mems+1));
        ERRTD2_timemem(e,t,z)=mean(TD2(:,e)-TD2(:,num_mems+1));  
        ERRWIND10_timemem(e,t,z)=mean(WIND10(:,e)-WIND10(:,num_mems+1));        
        
        MEANDIFFT2_timemem(e,t,z)=mean(modobsT2(:,1)-T2(:,e));       
        MEANDIFFTD2_timemem(e,t,z)=mean(modobsTD2(:,1)-TD2(:,e));
        MEANDIFFWIND10_timemem(e,t,z)=mean(modobsWIND10(:,1)-WIND10(:,e));
        
    end
    
    
    %%% RMSE mean, Variance %%%
    RMST2_all = sqrt((modobsT2(:,1)-modobsT2(:,2)).^2); RMST2(t,z) = mean(RMST2_all);
    biasT2(t,z) = mean(modobsT2(:,1)) - mean(modobsT2(:,2));
%     varT2(t)=mean(var(T2(:,1:num_mems)));
    
    RMSTD2_all = sqrt((modobsTD2(:,1)-modobsTD2(:,2)).^2); RMSTD2(t,z) = mean(RMSTD2_all);
    biasTD2(t,z) = mean(modobsTD2(:,1)) - mean(modobsTD2(:,2));
%     varTD2(t)=mean(var(TD2(:,1:num_mems)));
    
    RMSWIND10_all = sqrt((modobsWIND10(:,1)-modobsWIND10(:,2)).^2); RMSWIND10(t,z) = mean(RMSWIND10_all);   
    biasWIND10(t,z) = mean(modobsWIND10(:,1)) - mean(modobsWIND10(:,2));
%     varWIND10(t)=mean(var(WIND10(:,1:num_mems)));
    
    
    %%% Rank Histogram, Error to Variance, Variance
    for a=1:size(T2,1)
        obs_err_rand=normrnd(0,sigma_T2,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(T2(a,:)+obs_err_rand);
        rank=find(sorted==T2(a,num_mems+1));
        rank_sum_T2(rank,t)=rank_sum_T2(rank,t)+1;
        mod_var(a)=var(T2(a,1:num_mems));
    end
    obs_var(1:a)=sigma_T2^2;
    tot_varT2=obs_var+mod_var;
    errtovarT2(t,z)=( (modobsT2(:,2)-modobsT2(:,1))' * (modobsT2(:,2)-modobsT2(:,1))) / (sum(tot_varT2));
    varT2(t,z)=mean(mod_var);
    clearvars mod_var obs_var a
    
    
    for a=1:size(TD2,1)
        obs_err_rand=normrnd(0,sigma_TD2,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(TD2(a,:)+obs_err_rand);
        rank=find(sorted==TD2(a,num_mems+1));
        rank_sum_TD2(rank,t)=rank_sum_TD2(rank,t)+1;
        mod_var(a)=var(TD2(a,1:num_mems));
    end
    obs_var(1:a)=sigma_TD2^2;
    tot_varTD2=obs_var+mod_var;
    errtovarTD2(t,z)=( (modobsTD2(:,2)-modobsTD2(:,1))' * (modobsTD2(:,2)-modobsTD2(:,1))) / (sum(tot_varTD2));
    varTD2(t,z)=mean(mod_var);
    clearvars mod_var obs_var a
    
    
    for a=1:size(WIND10,1)
        obs_err_rand=normrnd(0,sigma_WIND10,[1,num_mems]); obs_err_rand(1,num_mems+1)=0;
        sorted=sort(WIND10(a,:)+obs_err_rand);
        rank=find(sorted==WIND10(a,num_mems+1));
        rank_sum_WIND10(rank,t)=rank_sum_WIND10(rank,t)+1;
        mod_var(a)=var(WIND10(a,1:num_mems));
    end
    obs_var(1:a)=sigma_WIND10^2;
    tot_varWIND10=obs_var+mod_var;
    errtovarWIND10(t,z)=( (modobsWIND10(:,2)-modobsWIND10(:,1))' * (modobsWIND10(:,2)-modobsWIND10(:,1))) / (sum(tot_varWIND10));
    varWIND10(t,z)=mean(mod_var);
    clearvars mod_var obs_var a
       
%     clearvars T2 TD2 U10 V10 WIND10 modobsT2 modobsTD2 modobsU10 modobsV10 modobsWIND10 
    %%% Build RMSE and total variance for all vars
    
    if strcmp(type,'fix')==1 && strcmp(d1fov,'fix')==1
        RMST2_all_sumcase(size(RMST2_all_sumcase,1)+1:size(RMST2_all_sumcase,1)+size(RMST2_all,1),1) = RMST2_all;
        tot_varT2_sumcase(size(tot_varT2_sumcase,1)+1:size(tot_varT2_sumcase,1)+size(tot_varT2,2),1)= tot_varT2;

        RMSTD2_all_sumcase(size(RMSTD2_all_sumcase,1)+1:size(RMSTD2_all_sumcase,1)+size(RMSTD2_all,1),1) = RMSTD2_all;
        tot_varTD2_sumcase(size(tot_varTD2_sumcase,1)+1:size(tot_varTD2_sumcase,1)+size(tot_varTD2,2),1)= tot_varTD2;

        RMSWIND10_all_sumcase(size(RMSWIND10_all_sumcase,1)+1:size(RMSWIND10_all_sumcase,1)+size(RMSWIND10_all,1),1) = RMSWIND10_all;
        tot_varWIND10_sumcase(size(tot_varWIND10_sumcase,1)+1:size(tot_varWIND10_sumcase,1)+size(tot_varWIND10,2),1)= tot_varWIND10;    
    end

    
end %cases




%%% Now that we've gone through each case for forecast hour (t)
%%% Apply bootstrap resampling to calculate sig bounds at that time %%%%%%%
if strcmp(type,'fix')==1 && strcmp(d1fov,'fix')==1
resample = 1000;

rs_etvT2 = zeros(1,resample);
samplelength = length(RMST2_all_sumcase);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errT2 = RMST2_all_sumcase(resample_inds).^2;   rs_varT2 = tot_varT2_sumcase(resample_inds);

rs_etvTD2 = zeros(1,resample);
samplelength = length(RMSTD2_all_sumcase);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errTD2 = RMSTD2_all_sumcase(resample_inds).^2;   rs_varTD2 = tot_varTD2_sumcase(resample_inds);

rs_etvWIND10 = zeros(1,resample);
samplelength = length(RMSWIND10_all_sumcase);
resample_inds = ceil(rand(samplelength,resample)*samplelength);
rs_errWIND10 = RMSWIND10_all_sumcase(resample_inds).^2;   rs_varWIND10 = tot_varWIND10_sumcase(resample_inds);


for i=1:resample;
    % height=[1000,925,850,700,500,400,300,250,200,150];
    rs_etvT2(i) = sum(rs_errT2(:,i)) / sum(rs_varT2(:,i));
    rs_etvTD2(i) = sum(rs_errTD2(:,i)) / sum(rs_varTD2(:,i));    
    rs_etvWIND10(i) = sum(rs_errWIND10(:,i)) / sum(rs_varWIND10(:,i));    
end

boot_etvT2 = sort(rs_etvT2);
boot_etvTD2 = sort(rs_etvTD2);
boot_etvWIND10 = sort(rs_etvWIND10);

siglevel=0.05;  ninetyith = siglevel*resample;
display(['Statistical Significance at bound: ' num2str(1-siglevel)])

sigbound_etvT2(t) = boot_etvT2(ninetyith);
sigbound_etvTD2(t) = boot_etvTD2(ninetyith);
sigbound_etvWIND10(t) = boot_etvWIND10(ninetyith);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



end %times (forecast hour)
end %domain


%%% Get Hourly correct mean of RMSE, ETV, VAR
for t=1:num_times
    %%% Temperature
    blahtotal=sum(obs_count_T2(t,:));
    obs_count_T2_rel(t,:)=obs_count_T2(t,:)/blahtotal;
 
    errtovarT2_time(t) = errtovarT2(t,:)*obs_count_T2_rel(t,:)';

    varT2_time(t) = varT2(t,:)*obs_count_T2_rel(t,:)';
    RMST2_time(t) = RMST2(t,:)*obs_count_T2_rel(t,:)';
    biasT2_time(t) = biasT2(t,:)*obs_count_T2_rel(t,:)';
    
    %%% Dewpoint Temperature
    blahtotal=sum(obs_count_TD2(t,:));
    obs_count_TD2_rel(t,:)=obs_count_TD2(t,:)/blahtotal;
    
    errtovarTD2_time(t) = errtovarTD2(t,:)*obs_count_TD2_rel(t,:)';
    varTD2_time(t) = varTD2(t,:)*obs_count_TD2_rel(t,:)';
    RMSTD2_time(t) = RMSTD2(t,:)*obs_count_TD2_rel(t,:)';
    biasTD2_time(t) = biasTD2(t,:)*obs_count_TD2_rel(t,:)';
    
    %%% Wind
    blahtotal=sum(obs_count_WIND10(t,:));
    obs_count_WIND10_rel(t,:)=obs_count_WIND10(t,:)/blahtotal;
    
    errtovarWIND10_time(t) = errtovarWIND10(t,:)*obs_count_WIND10_rel(t,:)';
    varWIND10_time(t) = varWIND10(t,:)*obs_count_WIND10_rel(t,:)';
    RMSWIND10_time(t) = RMSWIND10(t,:)*obs_count_WIND10_rel(t,:)';
    biasWIND10_time(t) = biasWIND10(t,:)*obs_count_WIND10_rel(t,:)';
    
    
    %%% Ensemble members
    RMST2_mem(:,t) = squeeze(RMST2_timemem(:,t,:))*obs_count_T2_rel(t,:)';
    RMSTD2_mem(:,t)= squeeze(RMSTD2_timemem(:,t,:))*obs_count_TD2_rel(t,:)';
    RMSWIND10_mem(:,t)=squeeze(RMSWIND10_timemem(:,t,:))*obs_count_WIND10_rel(t,:)';

    ERRT2_mem(:,t) = squeeze(ERRT2_timemem(:,t,:))*obs_count_T2_rel(t,:)';
    ERRTD2_mem(:,t)= squeeze(ERRTD2_timemem(:,t,:))*obs_count_TD2_rel(t,:)';
    ERRWIND10_mem(:,t)=squeeze(ERRWIND10_timemem(:,t,:))*obs_count_WIND10_rel(t,:)';    
  
    MEANDIFFT2_mem(:,t) = squeeze(MEANDIFFT2_timemem(:,t,:))*obs_count_T2_rel(t,:)';
    MEANDIFFTD2_mem(:,t) = squeeze(MEANDIFFTD2_timemem(:,t,:))*obs_count_TD2_rel(t,:)';
    MEANDIFFWIND10_mem(:,t) = squeeze(MEANDIFFWIND10_timemem(:,t,:))*obs_count_WIND10_rel(t,:)';
    
    
end

for z=1:size(list,1)
    blahtotal=sum(obs_count_T2(:,z));
    obs_count_T2_rel_case(:,z)=obs_count_T2(:,z)/blahtotal;
    blahtotal=sum(obs_count_TD2(:,z));
    obs_count_TD2_rel_case(:,z)=obs_count_TD2(:,z)/blahtotal;    
    blahtotal=sum(obs_count_WIND10(:,z));
    obs_count_WIND10_rel_case(:,z)=obs_count_WIND10(:,z)/blahtotal;  
    
    
    errtovarT2_case(z,1) = errtovarT2(:,z)'*obs_count_T2_rel_case(:,z);
    errtovarTD2_case(z,1) = errtovarTD2(:,z)'*obs_count_TD2_rel_case(:,z);
    errtovarWIND10_case(z,1) = errtovarWIND10(:,z)'*obs_count_WIND10_rel_case(:,z);
end

if strcmp(type,'fix')==1
fix_rank_sum_T2_accum=sum(rank_sum_T2,2);
fix_rank_sum_TD2_accum=sum(rank_sum_TD2,2);
fix_rank_sum_WIND10_accum=sum(rank_sum_WIND10,2);
elseif strcmp(type,'var')==1
    
var_rank_sum_T2_accum=sum(rank_sum_T2,2);
var_rank_sum_TD2_accum=sum(rank_sum_TD2,2);
var_rank_sum_WIND10_accum=sum(rank_sum_WIND10,2);
end
% diff_rank_sum_T2_accum=varied_rank_sum_T2_accum-fixed_rank_sum_T2_accum;
% diff_rank_sum_TD2_accum=varied_rank_sum_TD2_accum-fixed_rank_sum_TD2_accum;
% diff_rank_sum_WIND10_accum=varied_rank_sum_WIND10_accum-fixed_rank_sum_WIND10_accum;

 

figure
bar(eval([type '_rank_sum_T2_accum']))
xlim([1 num_mems+1])
title(['T2 Rank Histogram' domain ' ' type ' physics'])

figure
bar(eval([type '_rank_sum_TD2_accum']))
xlim([1 num_mems+1])
title(['TD2 Rank Histogram' domain ' ' type ' physics'])

figure
bar(eval([type '_rank_sum_WIND10_accum']))
xlim([1 num_mems+1])
title(['WIND10 Rank Histogram' domain ' ' type ' physics'])

cd(start)



set(groot,'defaultLineLineWidth',2)


%For domain 1 plots
%%%ETV
time=[0:48];
ideal=ones(49,1);
zero=zeros(49,1);
figure('paperunits','in','paperposition',[.5 .5 16 8])
% subplot(1,2,1)
h1=plot(time,D1fix_errtovarWIND10_time(:,1),'b'); hold on; plot(time,D1var_errtovarWIND10_time(:,1),'b--'); hold on; %plot(etvheight_D2WINDVAR2,height,'b:')
hold on; h4=plot(time,sigbound_etvWIND10,'LineWidth',1,'color','b');
hold on
% set(gca,'ydir','reverse')
h2=plot(time,D1fix_errtovarT2_time(:,1),'r'); hold on; plot(time,D1var_errtovarT2_time(:,1),'r--'); hold on; plot(time,sigbound_etvT2,'LineWidth',1,'color','r')
hold on
h3=plot(time,D1fix_errtovarTD2_time(:,1),'g'); hold on; plot(time,D1var_errtovarTD2_time(:,1),'g--'); hold on; plot(time,sigbound_etvTD2,'LineWidth',1,'color','g')
% plot(time,etvheight_D3WINDFIX,height,'r'); hold on; plot(etvheight_D3WINDVAR1,height,'r--'); hold on; plot(etvheight_D3WINDVAR2,height,'r:') 
hold on; plot(time,ideal,'k');set(gca,'FontSize',16)
hold on
h_legend=legend([h1 h2 h3 h4], {'10m Wind','2m T','2m TD', '95% sig'},'Location','NorthEast','fontsize', 20); hold on % h4, '95% sig'
set(h_legend, 'FontSize',16)
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics','fontsize', 20); hold on; axis tight; 
ylabel('Error to variance','fontsize', 20)
xlabel('Forecast hour','fontsize', 20)
xlim([0 48])
ylim([0.5 3])
hold on
set(gcf,'units','normalized','outerposition',[0 0 1 1])
% print('D1or2############_fix_var_sfc_etv_.png','-dpng','-r400')








%For domain 1 plots
%%%ETV
time=[0:48];
ideal=ones(49,1);
zero=zeros(49,1);
set(groot,'defaultLineLineWidth',2)
%%%RMSE
figure('paperunits','in','paperposition',[.5 .5 16 8])
subplot(3,1,1)
h1=plot(time,D1fix_RMSWIND10_time(:,1),'b'); hold on; plot(time,D1var_RMSWIND10_time(:,1),'b--');set(gca,'FontSize',16)% hold on; plot(etvheight_D2WINDVAR2,height,'b:')
hold on
h2=plot(time,D1fix_RMST2_time(:,1),'r'); hold on; plot(time,D1var_RMST2_time(:,1),'r--'); 
hold on
h3=plot(time,D1fix_RMSTD2_time(:,1),'g'); hold on; plot(time,D1var_RMSTD2_time(:,1),'g--'); 
hold on
legend([h1 h2 h3], {'10m Wind','2m T','2m TD'},'Location','NorthWest','fontsize', 12); hold on
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics'); hold on; axis tight; 
ylabel('RMSE','fontsize', 16)
xlabel('Forecast hour','fontsize', 16)
xlim([0 48])
ylim([1 3.5]) %3])

%%%variance
% figure('paperunits','in','paperposition',[.5 .5 16 8])
subplot(3,1,2)
h1=plot(time,D1fix_varWIND10_time(:,1),'b'); hold on; plot(time,D1var_varWIND10_time(:,1),'b--');set(gca,'FontSize',16)% hold on; plot(etvheight_D2WINDVAR2,height,'b:')
hold on
h2=plot(time,D1fix_varT2_time(:,1),'r'); hold on; plot(time,D1var_varT2_time(:,1),'r--'); 
hold on
h3=plot(time,D1fix_varTD2_time(:,1),'g'); hold on; plot(time,D1var_varTD2_time(:,1),'g--'); 
hold on
legend([h1 h2 h3], {'10m Wind','2m T','2m TD'},'Location','NorthWest','fontsize', 12); hold on
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics'); hold on; axis tight; 
ylabel('Variance','fontsize', 16)
xlabel('Forecast hour','fontsize', 16)
xlim([0 48])
ylim([0 6]) %5])

%%%bias
% figure('paperunits','in','paperposition',[.5 .5 16 8])
subplot(3,1,3)
h1=plot(time,D1fix_BIASWIND10_time(:,1),'b'); hold on; plot(time,D1var_BIASWIND10_time(:,1),'b--');set(gca,'FontSize',16)% hold on; plot(etvheight_D2WINDVAR2,height,'b:')
hold on
h2=plot(time,D1fix_BIAST2_time(:,1),'r'); hold on; plot(time,D1var_BIAST2_time(:,1),'r--'); 
hold on
h3=plot(time,D1fix_BIASTD2_time(:,1),'g'); hold on; plot(time,D1var_BIASTD2_time(:,1),'g--'); 
hold on
legend([h1 h2 h3], {'10m Wind','2m T','2m TD'},'Location','NorthWest','fontsize', 12); hold on
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics'); hold on; axis tight; 
hold on; plot(time,zero,'k')
ylabel('Bias','fontsize', 16)
xlabel('Forecast hour','fontsize', 16)
xlim([0 48])
ylim([-1 2.5]) %ylim([-0.5 2])

set(gcf,'units','normalized','outerposition',[0 0 1 1])
hold on
% annotation('textbox', [0 0.9 1 0.1],'FontSize',20,'String', 'Domain 1 Fixed (solid) vs Varied (dashed) physics','EdgeColor', 'none','HorizontalAlignment', 'center')
annotation('textbox', [0 0.9 1 0.1],'FontSize',20,'String', 'Domain 2 00 UTC (solid) vs 12 UTC (dashed)','EdgeColor', 'none','HorizontalAlignment', 'center')
% print('D1or2############_fix_var_sfc_allbutetv_.png','-dpng','-r400')







ddddddddddeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
set(groot,'defaultLineLineWidth',3)
time=[0:48];
ideal=ones(49,1);
figure('paperunits','in','paperposition',[.5 .5 16 8])
% subplot(1,2,1)
h1=plot(time,D1fixD2fix_errtovarWIND10_time(:,1),'b'); hold on; plot(time,D1varD2fix_errtovarWIND10_time(:,1),'b--'); hold on; plot(time,D1varD2var_errtovarWIND10_time(:,1),'b:');
h4=plot(time,sigbound_etvWIND10,'LineWidth',1,'color','b');
hold on
h2=plot(time,D1fixD2fix_errtovarT2_time(:,1),'r'); hold on; plot(time,D1varD2fix_errtovarT2_time(:,1),'r--'); hold on; plot(time,D1varD2var_errtovarT2_time(:,1),'r:'); hold on; plot(time,sigbound_etvT2,'LineWidth',1,'color','r');
hold on
h3=plot(time,D1fixD2fix_errtovarTD2_time(:,1),'g'); hold on; plot(time,D1varD2fix_errtovarTD2_time(:,1),'g--'); hold on; plot(time,D1varD2var_errtovarTD2_time(:,1),'g:'); hold on; plot(time,sigbound_etvTD2,'LineWidth',1,'color','g');
hold on; plot(time,ideal,'k')
hold on
legend([h1 h2 h3 h4], {'10m Wind','2m T','2m TD', '95% sig'},'Location','NorthWest','fontsize', 20); hold on
set(legend, 'FontSize',16)
% title('Domain 2 Fixed (D1 Fixed:solid; D1 Varied:dashed) vs Varied (D1 Varied:dotted) physics');
ylabel('Error to variance','fontsize', 20)
xlabel('Forecast hour','fontsize', 20)
xlim([0 48])
ylim([0.5 3.5])

set(gcf,'units','normalized','outerposition',[0 0 1 1])

% print('D2_fix_var_etv_sfc.png','-dpng','-r400')









%%% THIS WORKS FOR D1 FULL VAR, MP AND PBL ONLY VAR %%%

time=[0:48];
ideal=ones(49,1);
figure
h1=plot(time,D1var_errtovarWIND10_time(:,1),'b'); hold on; plot(time,D1PBL_errtovarWIND10_time(:,1),'b--'); hold on; plot(time,D1MP_errtovarWIND10_time(:,1),'b:');
hold on
h2=plot(time,D1var_errtovarT2_time(:,1),'r'); hold on; plot(time,D1PBL_errtovarT2_time(:,1),'r--'); hold on; plot(time,D1MP_errtovarT2_time(:,1),'r:');
hold on
h3=plot(time,D1var_errtovarTD2_time(:,1),'g'); hold on; plot(time,D1PBL_errtovarTD2_time(:,1),'g--'); hold on; plot(time,D1MP_errtovarTD2_time(:,1),'g:');
hold on; plot(time,ideal,'k')
hold on
legend([h1 h2 h3], {'10m Wind','2m T','2m TD'},'Location','NorthEast'); hold on
set(legend, 'FontSize',16)
title('Domain 1 varied physics(full varied:solid; MP only:dashed; PBL only:dotted)');
ylabel('Error-to-Variance')
xlabel('Forecast hour')
xlim([0 48])
% ylim([0 4])
% print('D2_fix_var_etv_sfc.png','-dpng','-r400')






%%%Compare on D2: varvar(full), MP and PBL var, and sigbound from fixfix?
%%%or fixfix itself?
time=[0:48];
set(groot,'defaultLineLineWidth',3)
ideal=ones(49,1);

figure('paperunits','in','paperposition',[.5 .5 16 8])
h1=plot(time,D1varD2var_errtovarWIND10_time(:,1),'b'); hold on; plot(time,D2PBL_errtovarWIND10_time(:,1),'b--'); hold on; plot(time,D2MP_errtovarWIND10_time(:,1),'b:'); hold on; plot(time,D1fixD2fix_errtovarWIND10_time(:,1),'LineWidth',1,'color','b');
hold on; set(gca,'FontSize',16); hold on
h2=plot(time,D1varD2var_errtovarT2_time(:,1),'r'); hold on; plot(time,D2PBL_errtovarT2_time(:,1),'r--'); hold on; plot(time,D2MP_errtovarT2_time(:,1),'r:'); hold on; plot(time,D1fixD2fix_errtovarT2_time(:,1),'LineWidth',1,'color','r');
hold on
h3=plot(time,D1varD2var_errtovarTD2_time(:,1),'g'); hold on; plot(time,D2PBL_errtovarTD2_time(:,1),'g--'); hold on; plot(time,D2MP_errtovarTD2_time(:,1),'g:');hold on; plot(time,D1fixD2fix_errtovarTD2_time(:,1),'LineWidth',1,'color','g');
hold on; plot(time,ideal,'k')
set(legend, 'FontSize',16)
legend([h1 h2 h3], {'10m Wind','2m T','2m TD'},'Location','Northwest','fontsize', 16); hold on
% title('Domain 2 Fixed (D1 Fixed:solid; D1 Varied:dashed) vs Varied (D1 Varied:dotted) physics');
ylabel('Error to variance','fontsize', 16)
xlabel('Forecast hour','fontsize', 16)
xlim([0 48])
ylim([0.5 2.5])
set(gcf,'units','normalized','outerposition',[0 0 1 1])




time=[0:48];
ideal=ones(49,1);
zero=zeros(49,1);
set(groot,'defaultLineLineWidth',2)
%%%RMSE
figure('paperunits','in','paperposition',[.5 .5 16 8])
subplot(3,1,1)
h1=plot(time,D1varD2var_RMSWIND10_time(:,1),'b'); hold on; plot(time,D2PBL_RMSWIND10_time(:,1),'b--'); hold on; plot(time,D2MP_RMSWIND10_time(:,1),'b:'); hold on; plot(time,D1fixD2fix_RMSWIND10_time(:,1),'LineWidth',1,'color','b');
hold on; set(gca, 'FontSize',14);
h2=plot(time,D1varD2var_RMST2_time(:,1),'r'); hold on; plot(time,D2PBL_RMST2_time(:,1),'r--'); hold on; plot(time,D2MP_RMST2_time(:,1),'r:'); hold on; plot(time,D1fixD2fix_RMST2_time(:,1),'LineWidth',1,'color','r');
hold on
h3=plot(time,D1varD2var_RMSTD2_time(:,1),'g'); hold on; plot(time,D2PBL_RMSTD2_time(:,1),'g--'); hold on; plot(time,D2MP_RMSTD2_time(:,1),'g:');hold on; plot(time,D1fixD2fix_RMSTD2_time(:,1),'LineWidth',1,'color','g');
hold on
legend([h1 h2 h3], {'10m Wind','2m T','2m TD'},'Location','NorthWest','fontsize', 20); hold on
set(legend, 'FontSize',12)
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics'); hold on; axis tight; 
ylabel('RMSE','fontsize', 16)
xlabel('Forecast hour','fontsize', 16)
xlim([0 48])
ylim([0.5 2.5]) %3])

%%%variance
% figure('paperunits','in','paperposition',[.5 .5 16 8])
subplot(3,1,2)
h1=plot(time,D1varD2var_varWIND10_time(:,1),'b'); hold on; plot(time,D2PBL_varWIND10_time(:,1),'b--'); hold on; plot(time,D2MP_varWIND10_time(:,1),'b:'); hold on; plot(time,D1fixD2fix_varWIND10_time(:,1),'LineWidth',1,'color','b');
hold on; set(gca, 'FontSize',14);
h2=plot(time,D1varD2var_varT2_time(:,1),'r'); hold on; plot(time,D2PBL_varT2_time(:,1),'r--'); hold on; plot(time,D2MP_varT2_time(:,1),'r:'); hold on; plot(time,D1fixD2fix_varT2_time(:,1),'LineWidth',1,'color','r');
hold on
h3=plot(time,D1varD2var_varTD2_time(:,1),'g'); hold on; plot(time,D2PBL_varTD2_time(:,1),'g--'); hold on; plot(time,D2MP_varTD2_time(:,1),'g:');hold on; plot(time,D1fixD2fix_varTD2_time(:,1),'LineWidth',1,'color','g');
hold on
% legend([h1 h2 h3], {'10m Wind','2m T','2m TD'},'Location','NorthWest'); hold on
set(legend, 'FontSize',12)
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics'); hold on; axis tight; 
ylabel('Variance','fontsize', 16)
xlabel('Forecast hour','fontsize', 16)
xlim([0 48])
ylim([0 6]) %5])

%%%bias
% figure('paperunits','in','paperposition',[.5 .5 16 8])
subplot(3,1,3)
h1=plot(time,D1varD2var_BIASWIND10_time(:,1),'b'); hold on; plot(time,D2PBL_BIASWIND10_time(:,1),'b--'); hold on; plot(time,D2MP_BIASWIND10_time(:,1),'b:'); hold on; plot(time,D1fixD2fix_BIASWIND10_time(:,1),'LineWidth',1,'color','b');
hold on; set(gca, 'FontSize',14);
h2=plot(time,D1varD2var_BIAST2_time(:,1),'r'); hold on; plot(time,D2PBL_BIAST2_time(:,1),'r--'); hold on; plot(time,D2MP_BIAST2_time(:,1),'r:'); hold on; plot(time,D1fixD2fix_BIAST2_time(:,1),'LineWidth',1,'color','r');
hold on
h3=blah;
plot(time,D1varD2var_BIASTD2_time(:,1),'g'); hold on; plot(time,D2PBL_BIASTD2_time(:,1),'g--'); hold on; plot(time,D2MP_BIASTD2_time(:,1),'g:');hold on; plot(time,D1fixD2fix_BIASTD2_time(:,1),'LineWidth',1,'color','g');
hold on
legend('fullvarvar','PBL','MP','fixfix')
% legend([h1 h2 h3], {'10m Wind','2m T','2m TD'},'Location','NorthWest','fontsize', 16); hold on
set(legend, 'FontSize',12)
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics'); hold on; axis tight; 
hold on; plot(time,zero,'k')
ylabel('Bias','fontsize', 16)
xlabel('Forecast hour','fontsize', 16)
xlim([0 48])
ylim([-0.5 1.5]) %ylim([-0.5 2])

set(gcf,'units','normalized','outerposition',[0 0 1 1])
hold on
% annotation('textbox', [0 0.9 1 0.1],'FontSize',20,'String', 'Domain 2 Fixed (D1 Fixed:solid; D1 Varied:dashed) vs Varied (D1 Varied:dotted) physics','EdgeColor', 'none','HorizontalAlignment', 'center')
% print('D1or2############_fix_var_sfc_allbutetv_.png','-dpng','-r400')
% print('D2_fix_var_sfc_allbutetv.png','-dpng','-r400')















%%%%%% ETV resulted over all forecast hours for each CASE %%%%%%%%%%%%%%
set(groot,'defaultLineLineWidth',2)
ideal=[1,1,1,1,1,1,1,1];
datelabels={'050600','050612','050700','050712','050800','050812','050900','050912'};
figure %figure('paperunits','in','paperposition',[.5 .5 16 8])
h1=plot(D1fixD2fix_etvWIND10,'b'); hold on; plot(D1varD2fix_etvWIND10,'b--'); hold on; plot(D1varD2var_etvWIND10,'b:');
hold on
h2=plot(D1fixD2fix_etvT2,'r'); hold on; plot(D1varD2fix_etvT2,'r--'); hold on; plot(D1varD2var_etvT2,'r:');
hold on
h3=plot(D1fixD2fix_etvTD2,'g'); hold on; plot(D1varD2fix_etvTD2,'g--'); hold on; plot(D1varD2var_etvTD2,'g:');
plot(ideal,'k'); hold on
hold on
legend([h1 h2 h3], {'10m Wind','2m T','2m TD'},'Location','NorthEast','fontsize', 16); hold on
set(legend, 'FontSize',16)
% title('Domain 2 Fixed (D1 Fixed:solid; D1 Varied:dashed) vs Varied (D1 Varied:dotted) physics');
ylabel('Error-to-Variance','fontsize', 16)
xlabel('Cases','fontsize', 16)
set(gca,'XtickLabel',datelabels,'XAxisLocation', 'bottom')
set(legend, 'FontSize',16)
set(gcf,'units','normalized','outerposition',[0 0 1 1])








time=[0:48];
ideal=ones(49,1);
figure
% subplot(1,2,1)
h1=plot(time,D1fixD2fix_varWIND10_time(:,1),'b'); hold on; plot(time,D1varD2fix_varWIND10_time(:,1),'b--'); hold on; plot(time,D1varD2var_varWIND10_time(:,1),'b:');
hold on
h2=plot(time,D1fixD2fix_varT2_time(:,1),'r'); hold on; plot(time,D1varD2fix_varT2_time(:,1),'r--'); hold on; plot(time,D1varD2var_varT2_time(:,1),'r:');
hold on
h3=plot(time,D1fixD2fix_varTD2_time(:,1),'g'); hold on; plot(time,D1varD2fix_varTD2_time(:,1),'g--'); hold on; plot(time,D1varD2var_varTD2_time(:,1),'g:');
hold on
legend([h1 h2 h3], {'10m Wind','2m T','2m TD'},'Location','NorthWest'); hold on
title('Domain 2 Fixed (D1 Fixed:solid; D1 Varied:dashed) vs Varied (D1 Varied:dotted) physics');
ylabel('Root Mean Square Error')
xlabel('Forecast hour')
xlim([0 48])
ylim([0 3])
% print('D2_fix_var_variance_sfc.png','-dpng','-r400')















%For domain 2 plots
time=[0:48];
ideal=ones(49,1);
zero=zeros(49,1);
set(groot,'defaultLineLineWidth',2)
%%%RMSE
figure('paperunits','in','paperposition',[.5 .5 16 8])
subplot(3,1,1)
h1=plot(time,D1fixD2fix_RMSWIND10_time(:,1),'b'); hold on; plot(time,D1varD2fix_RMSWIND10_time(:,1),'b--'); hold on; plot(time,D1varD2var_RMSWIND10_time(:,1),'b:')
hold on; set(legend, 'FontSize',14);
h2=plot(time,D1fixD2fix_RMST2_time(:,1),'r'); hold on; plot(time,D1varD2fix_RMST2_time(:,1),'r--'); hold on; plot(time,D1varD2var_RMST2_time(:,1),'r:');  
hold on
h3=plot(time,D1fixD2fix_RMSTD2_time(:,1),'g'); hold on; plot(time,D1varD2fix_RMSTD2_time(:,1),'g--'); hold on; plot(time,D1varD2var_RMSTD2_time(:,1),'g:'); 
hold on
legend([h1 h2 h3], {'10m Wind','2m T','2m TD'},'Location','NorthWest','fontsize', 20); hold on
set(legend, 'FontSize',12)
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics'); hold on; axis tight; 
ylabel('RMSE','fontsize', 16)
xlabel('Forecast hour','fontsize', 16)
xlim([0 48])
ylim([0.5 3]) %3])

%%%variance
% figure('paperunits','in','paperposition',[.5 .5 16 8])
subplot(3,1,2)
h1=plot(time,D1fixD2fix_varWIND10_time(:,1),'b'); hold on; plot(time,D1varD2fix_varWIND10_time(:,1),'b--'); hold on; plot(time,D1varD2var_varWIND10_time(:,1),'b:')
hold on; set(legend, 'FontSize',14);
h2=plot(time,D1fixD2fix_varT2_time(:,1),'r'); hold on; plot(time,D1varD2fix_varT2_time(:,1),'r--'); hold on; plot(time,D1varD2var_varT2_time(:,1),'r:');
hold on
h3=plot(time,D1fixD2fix_varTD2_time(:,1),'g'); hold on; plot(time,D1varD2fix_varTD2_time(:,1),'g--'); hold on; plot(time,D1varD2var_varTD2_time(:,1),'g:'); 
hold on
% legend([h1 h2 h3], {'10m Wind','2m T','2m TD'},'Location','NorthWest'); hold on
set(legend, 'FontSize',12)
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics'); hold on; axis tight; 
ylabel('Variance','fontsize', 16)
xlabel('Forecast hour','fontsize', 16)
xlim([0 48])
ylim([0 6]) %5])

%%%bias
% figure('paperunits','in','paperposition',[.5 .5 16 8])
subplot(3,1,3)
h1=plot(time,D1fixD2fix_BIASWIND10_time(:,1),'b'); hold on; plot(time,D1varD2fix_BIASWIND10_time(:,1),'b--'); hold on; plot(time,D1varD2var_BIASWIND10_time(:,1),'b:')
hold on; set(legend, 'FontSize',14);
h2=plot(time,D1fixD2fix_BIAST2_time(:,1),'r'); hold on; plot(time,D1varD2fix_BIAST2_time(:,1),'r--'); hold on; plot(time,D1varD2var_BIAST2_time(:,1),'r:');  
hold on
h3=plot(time,D1fixD2fix_BIASTD2_time(:,1),'g'); hold on; plot(time,D1varD2fix_BIASTD2_time(:,1),'g--'); hold on; plot(time,D1varD2var_BIASTD2_time(:,1),'g:');
hold on
% legend([h1 h2 h3], {'10m Wind','2m T','2m TD'},'Location','NorthWest','fontsize', 16); hold on
set(legend, 'FontSize',12)
% title('Domain 1 Fixed (solid) vs Varied (dashed) physics'); hold on; axis tight; 
hold on; plot(time,zero,'k')
ylabel('Bias','fontsize', 16)
xlabel('Forecast hour','fontsize', 16)
xlim([0 48])
ylim([-0.5 1.5]) %ylim([-0.5 2])

set(gcf,'units','normalized','outerposition',[0 0 1 1])
hold on
% annotation('textbox', [0 0.9 1 0.1],'FontSize',20,'String', 'Domain 2 Fixed (D1 Fixed:solid; D1 Varied:dashed) vs Varied (D1 Varied:dotted) physics','EdgeColor', 'none','HorizontalAlignment', 'center')
% print('D1or2############_fix_var_sfc_allbutetv_.png','-dpng','-r400')
% print('D2_fix_var_sfc_allbutetv.png','-dpng','-r400')



%%%%%%% For plotting observation locations on map projection %%%%%%%
        latlim=[22.0 52.0];
        lonlim=[-133.0 -77.0];
        hnd=figure;
        ax=worldmap(latlim,lonlim);
        latlim=getm(ax,'MapLatLimit');
        lonlim=getm(ax,'MapLonLimit');
        coast=load('coast');
        geoshow(ax, latlim, lonlim,'DisplayType', 'polygon', 'FaceColor', 'none')
        states=shaperead('usastatehi','UseGeocoords',true,'BoundingBox',[lonlim',latlim']);
        geoshow(ax,states,'FaceColor','none')
        close(hnd);
                
        scrsz = get(0,'ScreenSize');
        figure('Position',[1 scrsz(4)/1.5 scrsz(3)/1.5 scrsz(4)/1.5])
        hold on
        xlim(lonlim)
        ylim(latlim)
        plot(lonT2(in),latT2(in),'r+')
        plot(lonT2(~in),latT2(~in),'bo')

        geoshow(coast,'color','k')
        for is=1:size(states,1)
            plot(states(is).Lon,states(is).Lat,'k');    
        end
        
        


%%%%% Member plume plots for spread from ensemble and absolute error  %%%%%
time=[0:48];
ideal=ones(49,1);
zero=zeros(49,1);
set(groot,'defaultLineLineWidth',2)

time=[0:48];
zero=zeros(49,1);
figure
subplot(2,1,1)
h3 = plot(time,D1varD2var_MEANDIFFTD2_time_mem,'-','LineWidth',0.25); hold on; set(h3,'Color',[1 0.78 0.80]); hold on; h4 = plot(0,0,'Color',[1 0.78 0.80]); hold on
h2 = plot(time,D1varD2fix_MEANDIFFTD2_time_mem,'-','LineWidth',0.25); hold on; set(h2,'Color',[1 0.1 0]); hold on; h5 = plot(0,0,'Color',[1 0.1 0]); hold on
h1 = plot(time,D1fixD2fix_MEANDIFFTD2_time_mem,'-','LineWidth',0.25); hold on; set(h1,'Color',[0.2 0 0]); hold on; h6 = plot(0,0,'Color',[0.2 0 0]); hold on
hold on; plot(time,zero,'k','LineWidth',2); hold on;
set(gca,'FontSize',13)
xlabel('Forecast hour','fontsize',14)
ylabel('Member spread (m/s)','fontsize',14) %^{\circ}C
xlim([0 48])
% set(gca,'YTick',0:5:45);
ylim([-1.25 1.25])
set(gca,'YTick',-1.25:0.5:1.25);
legend([h6 h5 h4], {'D1FixD2Fix','D1VarD2Fix','D1VarD2Var'},'Location','NorthWest','fontsize',13);

subplot(2,1,2)

h3 = plot(time,D1varD2var_ERRTD2_time_mem,'-','LineWidth',0.25); hold on; set(h3,'Color',[1 0.78 0.80]); hold on; h4 = plot(0,0,'Color',[1 0.78 0.80]); hold on
h2 = plot(time,D1varD2fix_ERRTD2_time_mem,'-','LineWidth',0.25); hold on; set(h2,'Color',[1 0.1 0]); hold on; h5 = plot(0,0,'Color',[1 0.1 0]); hold on
h1 = plot(time,D1fixD2fix_ERRTD2_time_mem,'-','LineWidth',0.25); hold on; set(h1,'Color',[0.2 0 0]); hold on; h6 = plot(0,0,'Color',[0.2 0 0]); hold on
hold on; plot(time,zero,'k','LineWidth',2); hold on;
set(gca,'FontSize',13)
xlabel('Forecast hour','fontsize',14)
ylabel('Member error  (m/s)','fontsize',14)
xlim([0 48])
ylim([-1.5 1.5])
set(gca,'YTick',-1.5:0.5:1.5);
legend([h6 h5 h4], {'D1FixD2Fix','D1VarD2Fix','D1VarD2Var'},'Location','SouthWest','fontsize',13);





%%%compare pbl and mp only var with varvar (or fixfix?)
figure
subplot(2,1,1)
h3 = plot(time,D1varD2var_MEANDIFFTD2_time_mem,'-','LineWidth',0.25); hold on; set(h3,'Color',[1 0.78 0.80]); hold on; h4 = plot(0,0,'Color',[1 0.78 0.80]); hold on
h2 = plot(time,D2PBL_MEANDIFFTD2_time_mem,'-','LineWidth',0.25); hold on; set(h2,'Color',[1 0.1 0]); hold on; h5 = plot(0,0,'Color',[1 0.1 0]); hold on
h1 = plot(time,D2MP_MEANDIFFTD2_time_mem,'-','LineWidth',0.25); hold on; set(h1,'Color',[0.2 0 0]); hold on; h6 = plot(0,0,'Color',[0.2 0 0]); hold on
hold on; plot(time,zero,'k','LineWidth',2); hold on;
set(gca,'FontSize',13)
xlabel('Forecast hour','fontsize',14)
ylabel('Member spread (^{\circ}C)','fontsize',14) %^{\circ}C
xlim([0 48])
% set(gca,'YTick',0:5:45);
ylim([-1.25 1.25])
set(gca,'YTick',-1.25:0.5:1.25);
legend([h6 h5 h4], {'MP only','PBL only','Full var'},'Location','NorthWest','fontsize',13);

subplot(2,1,2)

h3 = plot(time,D1varD2var_ERRTD2_time_mem,'-','LineWidth',0.25); hold on; set(h3,'Color',[1 0.78 0.80]); hold on; h4 = plot(0,0,'Color',[1 0.78 0.80]); hold on
h2 = plot(time,D2PBL_ERRTD2_time_mem,'-','LineWidth',0.25); hold on; set(h2,'Color',[1 0.1 0]); hold on; h5 = plot(0,0,'Color',[1 0.1 0]); hold on
h1 = plot(time,D2MP_ERRTD2_time_mem,'-','LineWidth',0.25); hold on; set(h1,'Color',[0.2 0 0]); hold on; h6 = plot(0,0,'Color',[0.2 0 0]); hold on
hold on; plot(time,zero,'k','LineWidth',2); hold on;
set(gca,'FontSize',13)
xlabel('Forecast hour','fontsize',14)
ylabel('Member error  (^{\circ}C)','fontsize',14)
xlim([0 48])
ylim([-1.5 1.5])
set(gca,'YTick',-1:0.5:2.5);
legend([h6 h5 h4], {'MP only','PBL only','Full var'},'Location','SouthWest','fontsize',13);








%%%% Error of specific parameter option sub-ensembles 


%%%%%% Microphysics : 1-10=Thompson; 11-20=Morrison; 21-30=WSM6; 31-40=WDM6
thompson=[1:10];
morrison=[11:20];
wsm6=[21:30];
wdm6=[31:40];

time=[0:48];
ideal=ones(49,1);
zero=zeros(49,1);
set(groot,'defaultLineLineWidth',2)

figure('paperunits','in','paperposition',[.5 .5 16 8])
h1 = plot(time,D2PBL_ERRTD2_time_mem(thompson,:),'r','LineWidth',0.25); hold on; h5 = plot(0,0,'r'); hold on
h2 = plot(time,D2PBL_ERRTD2_time_mem(morrison,:),'g','LineWidth',0.25); hold on; h6 = plot(0,0,'g'); hold on
h3 = plot(time,D2PBL_ERRTD2_time_mem(wsm6,:),'b','LineWidth',0.25); hold on; h7 = plot(0,0,'b'); hold on
h4 = plot(time,D2PBL_ERRTD2_time_mem(wdm6,:),'m','LineWidth',0.25); hold on; h8 = plot(0,0,'m'); hold on

hold on; plot(time,zero,'k','LineWidth',2); hold on;
set(gca,'FontSize',16)
xlabel('Forecast hour','fontsize',16)
ylabel('Member error  (^{\circ}C)','fontsize',16) %^{\circ}C
xlim([0 48])
ylim([-.5 1.5])
% set(gca,'YTick',-1:0.5:2.5);
legend([h5 h6 h7 h8], {'Thompson','Morrison','WSM6','WDM6'},'Location','SouthWest','fontsize',25);
set(gcf,'units','normalized','outerposition',[0 0 1 1])




%%%%%%% PBL
ysu=[1,6,11,13,20,21,26,31,33,40];
myj=[2,5,12,14,19,22,25,32,34,39];
qnse=[3,8,9,15,18,23,28,29,35,38];
boul=[4,7,10,16,17,24,27,30,36,37];

time=[0:48];
ideal=ones(49,1);
zero=zeros(49,1);
set(groot,'defaultLineLineWidth',2)

figure('paperunits','in','paperposition',[.5 .5 16 8])
h1 = plot(time,D2PBL_ERRWIND10_time_mem(ysu,:),'r','LineWidth',0.25); hold on; h5 = plot(0,0,'r'); hold on
h2 = plot(time,D2PBL_ERRWIND10_time_mem(myj,:),'g','LineWidth',0.25); hold on; h6 = plot(0,0,'g'); hold on
h3 = plot(time,D2PBL_ERRWIND10_time_mem(qnse,:),'b','LineWidth',0.25); hold on; h7 = plot(0,0,'b'); hold on
h4 = plot(time,D2PBL_ERRWIND10_time_mem(boul,:),'m','LineWidth',0.25); hold on; h8 = plot(0,0,'m'); hold on

hold on; plot(time,zero,'k','LineWidth',2); hold on;
set(gca,'FontSize',16)
xlabel('Forecast hour','fontsize',16)
ylabel('Member error  (m/s)','fontsize',16) %^{\circ}C
xlim([0 48])
% ylim([-.5 1.5])
% set(gca,'YTick',-1:0.5:2.5);
legend([h5 h6 h7 h8], {'YSU','MYJ','QNSE','BouLac'},'Location','South','fontsize',25);
set(gcf,'units','normalized','outerposition',[0 0 1 1])






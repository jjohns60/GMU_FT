%% GMU-FT Application
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer: Jeremy Johnston
% Updated: April 3, 2021
%
% Description:
% This program facilitates the application of a soil freeze/thaw 
% classification model trained at George Mason University. The model is
% trained with thousands of globally distributed soil temperature sites, 
% remote sensing observations, and static land surface variables to acheive 
% global classification accuracies on the order of 88%. 
%
%
% Relevant publications:
%
% [1] Johnston, J., V. Maggioni, and P. Houser. 2019. “Investigating the 
%     Relationship Between Satellite-Based Freeze/Thaw Products and Land 
%     Surface Temperature.” IEEE Journal of Selected Topics in Applied 
%     Earth Observations and Remote Sensing, 1–25. 
%     https://doi.org/10.1109/JSTARS.2019.2926942.
%
% [2] Johnston, J., V. Maggioni, and P. Houser. 2020. “Comparing global 
%     passive microwave freeze/thaw records: Investigating differences 
%     between Ka- and L-band products.” Remote Sensing of Environment, 
%     247: 111936. https://doi.org/10.1016/j.rse.2020.111936.
%
% [3] Johnston, J. M., P. R. Houser, V. Maggioni, R. S. Kim, and C. 
%     Vuyovich. 2021. “Informing Improvements in Freeze/Thaw State 
%     Classification Using Subpixel Temperature.” IEEE Transactions on 
%     Geoscience and Remote Sensing, 1–19. 
%     https://doi.org/10.1109/TGRS.2021.3099292.
%
% [4] In production: Global Training and Performance
%
%-------------------------------------------------------------------------%
% Inputs:
% (1) Specify paths to data inputs for classification. Includes,
%
%   - SMAP Radiometer Twice-Daily rSIR-Enhanced EASE-Grid 2.0 Brightness 
%     Temperatures, Version 2
%  
%   - MODIS/Terra Snow Cover Daily L3 Global 0.05Deg CMG, Version 61
% 
%   - Rutgers Northern Hemisphere 24 km Weekly Snow Cover Extent, 
%     September 1980 Onward, Version 1
% 
% (2) Also, specify targeted start and end data of processing. This is in 
%     the form of a datetime array [dt_start dt_end] inclusive. If there is
%     not sufficient data to predict for a given data, will return outputs
%     as fill values (99)
%
%-------------------------------------------------------------------------%
% Outputs (2x daily, at 0.05 degree):
%   (1) Gridded binary FT classifications [frozen = 1, thawed = 0, 
%       fill = 99]
%   (2) Gridded probability of frozen [0% - 0 -> 100% - 1]
%   (3) Flags: 
%           0: no flag
%           1: probability of frozen between 40 - 60% (high uncertainty)
%           2: if not 1, probability of frozen between (25 - 75% (moderate
%              uncertainty)
%           3: ice cap (Koppen Climate Class 30), always frozen
%
%-------------------------------------------------------------------------%
%
% TO DO:
% (1) add functionality to variable resolution (currently 0.05-degree only)
% (2) improve initial thresholding
% (3) .....
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

%% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) input range of dates to model [datetime_start datetime_end]
processing_range = [datetime(2020,1,1) datetime(2020,12,31)];

%% (2) input paths
%save path
save_path = '/Volumes/GMU_FT/MODEL/OUTPUTS/';

%path to trained models
model_path = '/Volumes/GMU_FT/MODEL/TrainedModels/';

%path to folder containing all 8 ancillary inputs (.mat files at 0.05 grid)
static_predictors_path = '/Volumes/GMU_FT/MODEL/ancillary_files/';

%path to Rutgers snow product
GSL_path = '/Volumes/GMU_FT/MODEL/GSL_file/';

%path to MODIS NDSI data
NDSI_path = '/Volumes/GMU_FT/MODEL/NDSI_files/';

%path to land surface temperature data
LST_path = '/Volumes/GMU_FT/MODEL/LST_files/';

%path to SMAP 1.41 GHz brightness temperature data (3.125 km EASE grid)
SMAP_path = '/Volumes/GMU_FT/MODEL/SMAP_files/';

%path to SSMIS 19 GHz brightness temperature data (6.25 km EASE grid)
SSMIS19_path = '/Volumes/GMU_FT/MODEL/SSMIS_19files/';

%path to SSMIS 22 GHz brightness temperature data (6.25 km EASE grid)
SSMIS22_path = '/Volumes/GMU_FT/MODEL/SSMIS_22files/';

%path to SSMIS 37 GHz brightness temperature data (3.125 km EASE grid)
SSMIS37_path = '/Volumes/GMU_FT/MODEL/SSMIS_37files/';

%paths to relevant EASE2 grid information
EASE2_N3km_path = '/Volumes/GMU_FT/MODEL/grids/EASE2_N3.125km.geolocation.v0.9.nc';
EASE2_S3km_path = '/Volumes/GMU_FT/MODEL/grids/EASE2_S3.125km.geolocation.v0.9.nc';
EASE2_N6km_path = '/Volumes/GMU_FT/MODEL/grids/EASE2_N6.25km.geolocation.v0.9.nc';
EASE2_S6km_path = '/Volumes/GMU_FT/MODEL/grids/EASE2_S6.25km.geolocation.v0.9.nc';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%% MODEL APPLICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load in each regional FT model
model_files = dir([model_path 'RF*.mat']);
model_files = {model_files.name};
RF_Mdl_C10 = load([model_path model_files{contains(model_files,'Class10')}]);
n = fieldnames(RF_Mdl_C10); RF_Mdl_C10 = RF_Mdl_C10.(n{1});
RF_Mdl_C09 = load([model_path model_files{contains(model_files,'Class9')}]);
n = fieldnames(RF_Mdl_C09); RF_Mdl_C09 = RF_Mdl_C09.(n{1});
RF_Mdl_C08 = load([model_path model_files{contains(model_files,'Class8')}]);
n = fieldnames(RF_Mdl_C08); RF_Mdl_C08 = RF_Mdl_C08.(n{1});
RF_Mdl_C07 = load([model_path model_files{contains(model_files,'Class7')}]);
n = fieldnames(RF_Mdl_C07); RF_Mdl_C07 = RF_Mdl_C07.(n{1});
RF_Mdl_C06 = load([model_path model_files{contains(model_files,'Class6')}]);
n = fieldnames(RF_Mdl_C06); RF_Mdl_C06 = RF_Mdl_C06.(n{1});
RF_Mdl_C04 = load([model_path model_files{contains(model_files,'Class4')}]);
n = fieldnames(RF_Mdl_C04); RF_Mdl_C04 = RF_Mdl_C04.(n{1});
RF_Mdl_C03 = load([model_path model_files{contains(model_files,'Class3')}]);
n = fieldnames(RF_Mdl_C03); RF_Mdl_C03 = RF_Mdl_C03.(n{1});


%load in static predictors
staticVars = getStaticPredictors(static_predictors_path);

%Create data masks
%climate-based masking
climate_mask = staticVars.KC_Beck_10class == 1 | staticVars.KC_Beck_10class == 2 | ...
    staticVars.KC_Beck_10class == 5;
%masking of water dominated cells
water_mask = staticVars.water_proportion > 0.5;
%Always frozen mask
alwaysFZ_mask = staticVars.LC_MODIS_IGBP13 == 11;

%load in grid data to structures
EASE2_N3km.latitude = ncread(EASE2_N3km_path,'latitude'); EASE2_N3km.longitude = ncread(EASE2_N3km_path,'longitude');
EASE2_S3km.latitude = ncread(EASE2_S3km_path,'latitude'); EASE2_S3km.longitude = ncread(EASE2_S3km_path,'longitude');
EASE2_N6km.latitude = ncread(EASE2_N6km_path,'latitude'); EASE2_N6km.longitude = ncread(EASE2_N6km_path,'longitude');
EASE2_S6km.latitude = ncread(EASE2_S6km_path,'latitude'); EASE2_S6km.longitude = ncread(EASE2_S6km_path,'longitude');

%loop through each day for deriving 2x daily model outputs
dates_to_process = processing_range(1):processing_range(2);
for i = 220:length(dates_to_process)
    tic
    %get specific date
    date_i = dates_to_process(i);
    
    %identify files relevant to given date
    %(hemispheric (2) x overpass (2) x number of bands (n))
    SMAP_files = findFilesTB(SMAP_path,date_i); %SMAP 1.4GHz
    SSMIS19_files = findFilesTB(SSMIS19_path,date_i); %SSMIS 19GHz
    SSMIS22_files = findFilesTB(SSMIS22_path,date_i); %SSMIS 22GHz
    SSMIS37_files = findFilesTB(SSMIS37_path,date_i); %SSMIS 37GHz
    
    %extract data and re-grid to 0.05 degree resolution
    %6km products are fit to 0.1 degree grid, then resampled
    [M_22V,E_22V] = extractAndGrid(SSMIS22_path,SSMIS22_files,'-22V-',EASE2_N6km,EASE2_S6km,0.1); 
    M_22V = imresize(M_22V,2,'nearest'); E_22V = imresize(E_22V,2,'nearest');
    [M_19V,E_19V] = extractAndGrid(SSMIS19_path,SSMIS19_files,'-19V-',EASE2_N6km,EASE2_S6km,0.1); 
    M_19V = imresize(M_19V,2,'nearest'); E_19V = imresize(E_19V,2,'nearest');
    [M_19H,E_19H] = extractAndGrid(SSMIS19_path,SSMIS19_files,'-19H-',EASE2_N6km,EASE2_S6km,0.1); 
    M_19H = imresize(M_19H,2,'nearest'); E_19H = imresize(E_19H,2,'nearest');
    %3km products are fit to a 0.05 degree grid
    [M_37H,E_37H] = extractAndGrid(SSMIS37_path,SSMIS37_files,'-37H-',EASE2_N3km,EASE2_S3km,0.05); 
    [M_1H,E_1H] = extractAndGrid(SMAP_path,SMAP_files,'-1.4H-',EASE2_N3km,EASE2_S3km,0.05);
    [M_1V,E_1V] = extractAndGrid(SMAP_path,SMAP_files,'-1.4V-',EASE2_N3km,EASE2_S3km,0.05);
    
    %compute derived predictors
    %SMAP NPR
    M_NPR_SMAP = (M_1V - M_1H)./(M_1V + M_1H);
    E_NPR_SMAP = (E_1V - E_1H)./(E_1V + E_1H);
    %SSMIS NPR
    M_NPR_SSMIS = (M_19V - M_19H)./(M_19V + M_19H);
    E_NPR_SSMIS = (E_19V - E_19H)./(E_19V + E_19H);
    %37H - 19H spectral gradient
    M_SG_37H_19H = (M_37H - M_19H)./(37 - 19);
    E_SG_37H_19H = (E_37H - E_19H)./(37 - 19);
    %37H - 1.4H spectral gradient
    M_SG_37H_1H = (M_37H - M_1H)./(37 - 1.41);
    E_SG_37H_1H = (E_37H - E_1H)./(37 - 1.41);
    
    %get LST data
    [M_LST,E_LST] = loadLST(LST_path,date_i);

    %load in and re-grid snow cover data
    NH_SCE = getGSLSnow(GSL_path,date_i); %Northern hemisphere only (but gridded to full grid)
    NDSI = getNDSI(NDSI_path,date_i);
    
    
    %create table of all variables
    T = table();
    T.KC10 = staticVars.KC_Beck_10class(:);
    T.KC30_N_05 = staticVars.KC30_N_05(:);
    T.LC_MODIS_IGBP13 = staticVars.LC_MODIS_IGBP13(:);
    T.forest_proportion = staticVars.forest_proportion(:);
    T.water_proportion = staticVars.water_proportion(:);
    T.aspect = staticVars.aspect(:);
    T.TPI = staticVars.TPI(:);
    T.elev_sd_05 = staticVars.elev_sd_05(:);
    T.TB1_41V = M_1V(:);
    T.TB1_41H = M_1H(:);
    T.TB19V = M_19V(:);
    T.TB19H = M_19H(:);
    T.TB22V = M_22V(:);
    T.TB37H = M_37H(:);
    T.MODIS_LST = M_LST(:);
    T.MODIS_NDSI = NDSI(:);
    T.GSL_SCE = NH_SCE(:);
    T.NPR_SMAP = M_NPR_SMAP(:);
    T.NPR_SSMIS = M_NPR_SSMIS(:);
    T.SG_37H_19H = M_SG_37H_19H(:);
    T.SG_37H_1H = M_SG_37H_1H(:);
    
    %% predict with models (MORNRING-AM-OVERPASS) -> TO DO create single function later
    idx10 = staticVars.KC_Beck_10class == 10 & ~alwaysFZ_mask;
    [~,P10] = predict(RF_Mdl_C10,T(T.KC10 == 10 & ~alwaysFZ_mask(:),2:end));
    idx09 = staticVars.KC_Beck_10class == 9;
    [~,P09] = predict(RF_Mdl_C09,T(T.KC10 == 9,2:end));
    idx08 = staticVars.KC_Beck_10class == 8;
    [~,P08] = predict(RF_Mdl_C08,T(T.KC10 == 8,2:end));
    idx07 = staticVars.KC_Beck_10class == 7;
    [~,P07] = predict(RF_Mdl_C07,T(T.KC10 == 7,2:end));
    idx06 = staticVars.KC_Beck_10class == 6;
    [~,P06] = predict(RF_Mdl_C06,T(T.KC10 == 6,2:end));
    idx04 = staticVars.KC_Beck_10class == 4;
    [~,P04] = predict(RF_Mdl_C04,T(T.KC10 == 4,2:end));
    idx03 = staticVars.KC_Beck_10class == 3;
    [~,P03] = predict(RF_Mdl_C03,T(T.KC10 == 3,2:end));
    
    %combine into grid of all binart FT values
    M_FZp = NaN([3600 7200]);
    M_FZp(idx10) = P10(:,2);
    M_FZp(idx09) = P09(:,2);
    M_FZp(idx08) = P08(:,2);
    M_FZp(idx07) = P07(:,2);
    M_FZp(idx06) = P06(:,2);
    M_FZp(idx04) = P04(:,2);
    M_FZp(idx03) = P03(:,2);
    M_FZp(alwaysFZ_mask) = 1;
    M_FZp(water_mask) = NaN;
    
    %convert to binary FT state
    M_FT = NaN([3600 7200]);
    M_FT(M_FZp > 0.5) = 1;
    M_FT(M_FZp <= 0.5) = 0;
    
    %compute flags
    M_Flags = zeros([3600 7200]);
    M_Flags(M_FZp > 0.25 & M_FZp < 0.75) = 2; %flag for regions with FZp between 25 and 75%
    M_Flags(M_FZp > 0.4 & M_FZp < 0.6) = 1; %flag for regions with FZp between 40 and 60%
    M_Flags(alwaysFZ_mask) = 3; %flag for regions deemed always frozen
    M_Flags(climate_mask) = 4; %flag for regions deemed always thawed/fill (model is not applied)
    M_Flags(water_mask) = 5; %flag for water covered regions
    
    %Update table for afternoon overpass computations
    T.TB1_41V = E_1V(:);
    T.TB1_41H = E_1H(:);
    T.TB19V = E_19V(:);
    T.TB19H = E_19H(:);
    T.TB22V = E_22V(:);
    T.TB37H = E_37H(:);
    T.MODIS_LST = E_LST(:);
    T.NPR_SMAP = E_NPR_SMAP(:);
    T.NPR_SSMIS = E_NPR_SSMIS(:);
    T.SG_37H_19H = E_SG_37H_19H(:);
    T.SG_37H_1H = E_SG_37H_1H(:);
    
    %predict with models (AFTERNOON/EVENING-PM-OVERPASS)
    idx10 = staticVars.KC_Beck_10class == 10 & ~alwaysFZ_mask;
    [~,P10] = predict(RF_Mdl_C10,T(T.KC10 == 10 & ~alwaysFZ_mask(:),2:end));
    idx09 = staticVars.KC_Beck_10class == 9;
    [~,P09] = predict(RF_Mdl_C09,T(T.KC10 == 9,2:end));
    idx08 = staticVars.KC_Beck_10class == 8;
    [~,P08] = predict(RF_Mdl_C08,T(T.KC10 == 8,2:end));
    idx07 = staticVars.KC_Beck_10class == 7;
    [~,P07] = predict(RF_Mdl_C07,T(T.KC10 == 7,2:end));
    idx06 = staticVars.KC_Beck_10class == 6;
    [~,P06] = predict(RF_Mdl_C06,T(T.KC10 == 6,2:end));
    idx04 = staticVars.KC_Beck_10class == 4;
    [~,P04] = predict(RF_Mdl_C04,T(T.KC10 == 4,2:end));
    idx03 = staticVars.KC_Beck_10class == 3;
    [~,P03] = predict(RF_Mdl_C03,T(T.KC10 == 3,2:end));
    
    
    %combine into grid of all binary FT values
    E_FZp = NaN([3600 7200]);
    E_FZp(idx10) = P10(:,2);
    E_FZp(idx09) = P09(:,2);
    E_FZp(idx08) = P08(:,2);
    E_FZp(idx07) = P07(:,2);
    E_FZp(idx06) = P06(:,2);
    E_FZp(idx04) = P04(:,2);
    E_FZp(idx03) = P03(:,2);
    E_FZp(alwaysFZ_mask) = 1;
    E_FZp(water_mask) = NaN;
    
    %convert to binary FT state
    E_FT = NaN([3600 7200]);
    E_FT(E_FZp > 0.5) = 1;
    E_FT(E_FZp <= 0.5) = 0;
    
    %compute flags
    E_Flags = zeros([3600 7200]);
    E_Flags(E_FZp > 0.25 & E_FZp < 0.75) = 2; %flag for regions with FZp between 25 and 75%
    E_Flags(E_FZp > 0.4 & E_FZp < 0.6) = 1; %flag for regions with FZp between 40 and 60%
    E_Flags(alwaysFZ_mask) = 3; %flag for regions deemed always frozen
    E_Flags(climate_mask) = 4; %flag for regions deemed always thawed/fill (model is not applied)
    E_Flags(water_mask) = 5; %flag for water covered regions

    %save daily output to save path
    save([save_path 'GMU_FT_' num2str(year(date_i)) num2str(day(date_i,'dayofyear'),'%03d') '_0_05_V1.mat'],'M_FT','M_FZp','M_Flags','E_FT','E_FZp','E_Flags');

    disp([datestr(date_i) ' completed processing in ' num2str(toc) 'seconds'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












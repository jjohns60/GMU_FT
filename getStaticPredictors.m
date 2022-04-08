function S = getStaticPredictors(path)
%getStaticPredictors Takes an input path and loads in all relevant
%variables to a structure. Works if files are specifically named (see
%below)

%load in static variables and append to structure
static_files = dir([path '*.mat']);
idx = startsWith({static_files.name},'.');  %exclude hidden files
static_files = static_files(~idx);

%climate classes
KC_Beck_10class = load([path static_files(contains({static_files.name},'climate_10class')).name]);
names = fieldnames(KC_Beck_10class);
S.KC_Beck_10class = KC_Beck_10class.(names{1});

%Number of climate classes
KC30_N_05 = load([path static_files(contains({static_files.name},'Nclasses')).name]);
names = fieldnames(KC30_N_05);
S.KC30_N_05 = KC30_N_05.(names{1});

% land cover
LC_MODIS_IGBP13 = load([path static_files(contains({static_files.name},'landcover')).name]);
names = fieldnames(LC_MODIS_IGBP13);
S.LC_MODIS_IGBP13 = LC_MODIS_IGBP13.(names{1});

% forest proportion
forest_proportion = load([path static_files(contains({static_files.name},'forest')).name]);
names = fieldnames(forest_proportion);
S.forest_proportion = forest_proportion.(names{1});

% surface water proportion
water_proportion = load([path static_files(contains({static_files.name},'water')).name]);
names = fieldnames(water_proportion);
S.water_proportion = water_proportion.(names{1});

%aspect
aspect = load([path static_files(contains({static_files.name},'aspect')).name]);
names = fieldnames(aspect);
S.aspect = aspect.(names{1});

%TPI
TPI = load([path static_files(contains({static_files.name},'TPI')).name]);
names = fieldnames(TPI);
S.TPI = TPI.(names{1});

%elevation standard deviation
elev_sd_05 = load([path static_files(contains({static_files.name},'elevation_std')).name]);
names = fieldnames(elev_sd_05);
S.elev_sd_05 = elev_sd_05.(names{1});

end


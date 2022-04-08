function [M,E] = extractAndGrid(path,files,spec,Ngrid,Sgrid,res)
%extractAndGrid Takes an array of input files located on a given path, and
%identifies ones to process using a unique string/band (spec, i.e., -22V-).
%To speed processing, N/Sgrid is a structure containing information 
%S.latitude (latitudes) and S.longitude (longitudes). These correspond to 
%the data located in the input 'files'. Resolution of  re-gridding also 
%input as 'res' in degrees. Specified as both north specific (N) and south
%specific (S) polar grids
%
%
%derived from need to combine hemispheric dataset, therefor this program
%first identifies northern hemisphere data with the given overpass
%indicator (E,M), then southern. Outputs 2 matched resolution grids, one 
%for the morning (M, AM) period, and one for the afternoon period (E, PM)

%{
%testing
path = SSMIS22_path;
files = SSMIS22_files;
spec = '-22V-';
Ngrid = EASE2_N6km;
Sgrid = EASE2_S6km;
res = 0.1;
%}

%create normally spaces global grid using res
lat_target = (90-(0.5*res):-res:-90+(0.5*res))';
lon_target = (-180+(0.5*res):res:180-(0.5*res));
lat_target = repmat(lat_target,[1 length(lon_target)]);
lon_target = repmat(lon_target,[size(lat_target,1) 1]);

%identify relevant files by hemisphere and overpass time
Nfile_M = files(contains(files,'EASE2_N') & contains(files,[spec 'M']));
Sfile_M = files(contains(files,'EASE2_S') & contains(files,[spec 'M']));
Nfile_E = files(contains(files,'EASE2_N') & contains(files,[spec 'E']));
Sfile_E = files(contains(files,'EASE2_S') & contains(files,[spec 'E']));

%{
%get file information for formatting (not needed)
I = ncinfo([path Nfile_M{1}]); %assumes for same band/product that data will be formatted the same
fill_value = I.Variables(5).Attributes(4); fill_value = fill_value.Value;
missing_value = I.Variables(5).Attributes(5); missing_value = missing_value.Value;
valid_range = I.Variables(5).Attributes(6); valid_range = valid_range.Value;
%}


%load in band data for each (already in correct units)
%if no file exits, fill placeholder with NaNs
if length(Nfile_M) == 1
    %load in data
    Ndata_M = single(ncread([path Nfile_M{1}],'TB'));
    %regrid evening overpass data
    N_M = gridvaluesearch(Ndata_M,Ngrid.latitude,Ngrid.longitude,lat_target,lon_target,'average');
    %create mask
    m = isnan(N_M);
    %m(0.1*size(m,1),:) = 0;
    m = bwareaopen(m,0.5*(length(m(:))));
    %interpolate missing data horizontally
    N_M = fillmissing(N_M,'linear',2);
    N_M(m) = NaN;
else
    N_M = NaN(size(lat_target));
end

if length(Sfile_M) == 1
    %load in data
    Sdata_M = single(ncread([path Sfile_M{1}],'TB'));
    %regrid evening overpass data
    S_M = gridvaluesearch(Sdata_M,Sgrid.latitude,Sgrid.longitude,lat_target,lon_target,'average');
    %create mask
    m = isnan(S_M);
    %m(end-0.1*size(m,1),:) = 0;
    m = bwareaopen(m,0.5*(length(m(:))));
    %interpolate missing data horizontally
    S_M = fillmissing(S_M,'linear',2);
    S_M(m) = NaN;
else
    S_M = NaN(size(lat_target));
end

if length(Nfile_E) == 1
    %load in data
    Ndata_E = single(ncread([path Nfile_E{1}],'TB'));
    %regrid evening overpass data
    N_E = gridvaluesearch(Ndata_E,Ngrid.latitude,Ngrid.longitude,lat_target,lon_target,'average');
    %create mask
    m = isnan(N_E);
    %m(0.1*size(m,1),:) = 0;
    m = bwareaopen(m,0.5*(length(m(:))));
    %interpolate missing data horizontally
    N_E = fillmissing(N_E,'linear',2);
    N_E(m) = NaN;
else
    N_E = NaN(size(lat_target));
end

if length(Sfile_E) == 1
    %load in data
    Sdata_E = single(ncread([path Sfile_E{1}],'TB'));
    %regrid evening overpass data
    S_E = gridvaluesearch(Sdata_E,Sgrid.latitude,Sgrid.longitude,lat_target,lon_target,'average');
    %create mask
    m = isnan(S_E);
    %m(0.1*size(m,1),:) = 0;
    m = bwareaopen(m,0.5*(length(m(:))));
    %interpolate missing data horizontally
    S_E = fillmissing(S_E,'linear',2);
    S_E(m) = NaN;
else
    S_E = NaN(size(lat_target));
end


%recombine into a single global grid for morning (M)
M = N_M;
M(0.5*size(M,1)+1:end,:) = S_M(0.5*size(M,1)+1:end,:);

%and evening (E)
E = N_E;
E(0.5*size(E,1)+1:end,:) = S_E(0.5*size(E,1)+1:end,:);

end


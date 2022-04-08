function NH_SCE = getGSLSnow(path,date)
%getGSLSnow Get snow cover extent information from Rutgers GSL weekly
%product and re-grid to 0.05 degrees

%{
path = GSL_path;
date = datetime(2020,1,1);
%}


%identify file
file = dir([path 'G10035*.nc']);

%load data
d = ncread([path file.name],'time');
t = datetime(1980,8,25) + days(d);

%ensure date is valid
if date < t(1) || date > t(end)
    NH_SCE = NaN([3600 7200]);
    return
end
%identify snow cover layer to extract
dur = abs(days(date - t));
I = find(dur == min(dur));

%get snow cover data and re-grid
y = ncread([path file.name],'latitude');
x = ncread([path file.name],'longitude');
d = ncread([path file.name],'snow_cover_extent');
d = d(:,:,I); %extract relevant time step
d(d < 0) = NaN;

%gridding
res = 0.25;
lat_target = (90-(0.5*res):-res:-90+(0.5*res))';
lon_target = (-180+(0.5*res):res:180-(0.5*res));
lat_target = repmat(lat_target,[1 length(lon_target)]);
lon_target = repmat(lon_target,[size(lat_target,1) 1]);
NH_SCE = gridvaluesearch(d,y,x,lat_target,lon_target,'average');
%create mask
m = isnan(NH_SCE);
m = bwareaopen(m,0.5*(length(m(:))));
%interpolate missing data horizontally
NH_SCE = fillmissing(NH_SCE,'linear',2);
NH_SCE(m) = NaN;

%resample to 0.05 grid (from quarter degree)
NH_SCE = imresize(NH_SCE,5,'nearest');
NH_SCE(NH_SCE > 1) = 1;
NH_SCE(NH_SCE < 0) = NaN;

end


function [M,E] = loadLST(path,date)
%loadLST Extracts MODIS land surface temperature 0.05 data from Terra/Aqua
%products and combines this information into a single morning (M) and
%evening (E) overpass

%load in all filenames
files = dir([path 'M*.hdf']);
filenames = {files.name};

%create date string to find
date_str = ['A' num2str(year(date),'%04d') num2str(day(date,'dayofyear'),'%03d') '.'];

%identify and store names of all files containing date
filenames = filenames(contains(filenames,date_str));

if length(filenames) == 0
    M = NaN([3600 7200]);
    E = NaN([3600 7200]);
    return
end

try
%load in and format LST data (remove fill and use scale factor
LST_Day1 = double(hdfread([path filenames{1}],'LST_Day_CMG'));
LST_Day1(LST_Day1 == 0) = NaN;
LST_Day1 = LST_Day1 .* 0.02;
LST_Night1 = double(hdfread([path filenames{1}],'LST_Night_CMG'));
LST_Night1(LST_Night1 == 0) = NaN;
LST_Night1 = LST_Night1 .* 0.02;
catch
    LST_Night1 = NaN([3600 7200]);
    LST_Day1 = NaN([3600 7200]);
end

try
%second day
LST_Day2 = double(hdfread([path filenames{2}],'LST_Day_CMG'));
LST_Day2(LST_Day2 == 0) = NaN;
LST_Day2 = LST_Day2 .* 0.02;
LST_Night2 = double(hdfread([path filenames{2}],'LST_Night_CMG'));
LST_Night2(LST_Night2 == 0) = NaN;
LST_Night2 = LST_Night2 .* 0.02;
catch
    LST_Night2 = NaN([3600 7200]);
    LST_Day2 = NaN([3600 7200]);
end

%average where both values are present, otherwise retain one or the other
%(intended to maximize LST daily coverage/minimize loss due to cloud cover)
idx_vv = ~isnan(LST_Day1) & ~isnan(LST_Day2);
M = NaN(size(LST_Day1));
M(idx_vv) = (LST_Day1(idx_vv) + LST_Day2(idx_vv))./2;
idx_nv = isnan(LST_Day1) & ~isnan(LST_Day2);
M(idx_nv) = LST_Day2(idx_nv);
idx_vn = ~isnan(LST_Day1) & isnan(LST_Day2);
M(idx_vn) = LST_Day1(idx_vn);

idx_vv = ~isnan(LST_Night1) & ~isnan(LST_Night2);
E = NaN(size(LST_Night1));
E(idx_vv) = (LST_Night1(idx_vv) + LST_Night2(idx_vv))./2;
idx_nv = isnan(LST_Night1) & ~isnan(LST_Night2);
E(idx_nv) = LST_Night2(idx_nv);
idx_vn = ~isnan(LST_Night1) & isnan(LST_Night2);
E(idx_vn) = LST_Night1(idx_vn);

end
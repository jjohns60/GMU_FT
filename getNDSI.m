function NDSI = getNDSI(path,date)
%getNDSI Get and regrid relevant MODIS NDSI data within a specified path
%and for a specified date

%{
%testing
path = NDSI_path;
date = datetime(2020,1,15);
%}

%get all filenames
files = dir([path 'MYD*.hdf']);
filenames = {files.name};

%create date string to find
date_str = ['A' num2str(year(date),'%04d') num2str(day(date,'dayofyear'),'%03d') '.'];

%identify and store names of all files containing date
filenames = filenames(contains(filenames,date_str));

%extract data (already at 0.05 resolution)
if length(filenames) == 1
    NDSI = double(hdfread([path filenames{1}],'Day_CMG_Snow_Cover'));
    NDSI(NDSI < 0) = NaN;
    NDSI(NDSI > 100) = NaN;
    NDSI = NDSI./100;
else
    NDSI = NaN([3600 7200]);
end


end


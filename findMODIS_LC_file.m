function filepath = findMODIS_LC_file(path,date)
%findMODIS_LC_file Takes an input path to a folder containing MODIS land
%cover data (MCD12C1, 0.05 degree) and a date, and determines the file with 
%the year nearest to that of the input date

%{
%testing
path = '/Volumes/GMU_FT/DATA/ANCILLARY/MODIS_LANDCOVER/';
date = datetime(2021,1,1);
%}

%get year from date
current_year = year(date);

%get all file information (.hdf format only)
files = dir([path '*.hdf']);
files = files(~startsWith({files.name},'.'));
files = {files.name};

%get years associated with each file
file_years = NaN(size(files));
for i = 1:length(files)
    file = files{i};
    file_years(i) = str2double(file(10:13));
end

%get index of nearest filename
[~,idx] = min(abs(current_year - file_years));

%get and return the filename with attached path
filepath = [path files{idx}];

end


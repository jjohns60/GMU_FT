function filenames = findFilesTB(path,date)
%findFiles Searches an input path for all data from a given day.
%Specifically works for NSIDC EASE2 gridded TB files

%get all file names
files = dir([path 'NSIDC*.nc']);
filenames = {files.name};

%create date string to find
date_str = ['-' num2str(year(date),'%04d') num2str(day(date,'dayofyear'),'%03d') '-'];

%identify and store names of all files containing date
filenames = filenames(contains(filenames,date_str));
end


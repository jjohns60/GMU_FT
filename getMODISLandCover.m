function [LC_13,Pforest,Pwater] = getMODISLandCover(file)
%getMODISLandCover Takes an input MODIS file (as path) and returns a 13 
%class land cover, forested proportion, and water proportion at the file
%resolution for the input MODIS land cover (~0.05 degree grid)

%{
%test
file = '/Volumes/GMU_FT/DATA/ANCILLARY/MODIS_LANDCOVER/MCD12C1.A2020001.006.2021362215328.hdf';
%}

%load in data
LC = double(hdfread(file,'Majority_Land_Cover_Type_1'));
LC_pct = hdfread(file,'Land_Cover_Type_1_Percent');

%get proportion of water in each cell
Pwater = double(LC_pct(:,:,1));
Pwater(Pwater == 255) = NaN;
Pwater = Pwater./100;

%get estimated proportion of tree cover in each cell (estimated from class proportions)
Pforest = double(LC_pct(:,:,[2:6,9,10,13,15])); %extract tree relevant classes
Pforest(Pforest == 255) = NaN;
%weight each layer in grid by estimated proportion of trees
Pforest(:,:,1:5) = Pforest(:,:,1:5).*0.8; %forested classes
Pforest(:,:,6) = Pforest(:,:,6).*0.45; %woody savannas
Pforest(:,:,7) = Pforest(:,:,7).*0.2; %savannas
Pforest(:,:,8) = Pforest(:,:,8).*0.1; %croplands
Pforest(:,:,9) = Pforest(:,:,9).*0.25; %crop and natural veg mosaic
%normalize, 80% is maximum possible
Pforest = sum(Pforest,3)./80;

%convert to 13 classes
LC_13 = NaN(size(LC));
LC_13(LC == 0) = 13; %water = 13
LC_13(LC == 16) = 12; %barren = 12
LC_13(LC == 15) = 11; %snow and ice = 11
LC_13(LC == 13) = 10; %urban and built up = 10
LC_13(LC == 14 | LC == 12) = 9; %croplands = 9
LC_13(LC == 11) = 8; %wetlands = 8
LC_13(LC == 10) = 7; %grasslands = 7
LC_13(LC == 9 | LC == 8) = 6; %savannas = 6
LC_13(LC == 7 | LC == 6) = 5; %shrublands = 5
LC_13(LC == 5) = 4; %mixed forest = 4
LC_13(LC == 4 | LC == 3) = 3; %deciduous = 3
LC_13(LC == 2) = 2; %evergreen broadleaf = 2
LC_13(LC == 1) = 1; %evergreen needleleaf = 1


end


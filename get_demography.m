clc
clear all

%Import the location, area and total population of all neighborhoods in
%Jalisco state
raw_table = readtable('area_population.csv');

jalisco.pobtot = table2array(raw_table(:,3));
jalisco.area = table2array(raw_table(:,9));
jalisco.y = table2array(raw_table(:,10));
jalisco.x = table2array(raw_table(:,11));

%Calculate population density
jalisco.denspob = jalisco.pobtot ./ ( jalisco.area*1e-6);

%Define a window around the Guadalajara Metropolitan Area
GDL.xwindow = [6.5, 6.9]*1e5;
GDL.ywindow = [2.26 2.305]*1e6;

nancheck = ( jalisco.x > GDL.xwindow(1)) .* (jalisco.x < GDL.xwindow(2));
nancheck = double(nancheck).* ((jalisco.y>GDL.ywindow(1)).*(jalisco.y < GDL.ywindow(2)));
nancheck = logical(nancheck);

% Import census data showing population by age and healthcare coverage
cens2010 = readtable('inegi_jalisco.csv');

%Get lat lon as decimals
 cens.lat = [];

cens.lon = char(table2array(cens2010(:,4)));
cens.lon = str2num(cens.lon(:,1:3)) + ( str2num(cens.lon(:,4:5)) + ...
    str2num(cens.lon(:,6:7))/60)/60;


%for element = 1:size(cens2010,1)

 %   lon = cell2mat(table2array
    
    
%end

%Define grid to represent the city
load('malla892nod.mat');

%Calculate center of all grid elements
grid.centroid = [];

for element = 1:length(grid.mm)
    
    x = mean( grid.z(1,grid.mm(:,element)));
    y = mean( grid.z(2,grid.mm(:,element)));
    
    grid.centroid(element,:) = [x,y]; 
    
end







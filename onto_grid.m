clc
clear all

%Obtain demography for Guadalajara and map onto our numerical grid. We
%begin with raw variables from inegi census data products (2010) and our
%numerical grid. 


%Send your questions to:
% Noel Brizuela | nogutier@ucsd.edu
% Scripps Institution of Oceanography, University of California, San Diego
% December of 2018

%% LOAD RAW DATA WITH CITY BLOCK RESOLUTION
census = readtable('jalisco_2010_ageb_data.csv');

%Keep AGEB ID (1,3,5,7), total population (9), population 15 or older (24)
%and population with no health insurance (145)
census = census(:, [1,3,5,7,8,9,24,145]);

%Extract data aggregated by AGEB
ghostcheck = (table2array(census(:,5)) == 0);

census = census(ghostcheck,:);

ageb = char(census.ageb);

ghostcheck = NaN(size(census,1),1);

for block = 1:size(census,1)
    
    if sum(ageb(block,:) == char('0000')) == 4
        ghostcheck(block) = 0;
    else
        ghostcheck(block) = 1;
    end
        
end

%Data aggregated by AGEB we are capturing 6.399 M out of the 7.350 M people 
%who live in Jalisco state. The rest live outside urban AGEBs (i.e. in
%rural areas)
census = census(logical(ghostcheck),:);

clear ghostcheck block ageb

AGEB = ["h0"];

for k = 1:size(census,1)
    
    AGEB(k) = join(['14', sprintf('%03d',census.mun(k)), sprintf('%04d', ...
        census.loc(k)),  string(census.ageb(k))],'');
    
end

AGEB = AGEB';

census = [table(AGEB), census(:,6:end)];

clear AGEB k

save('AGEB_aggregated.mat');

%% LOAD AGGREGATED UNIT DATA

load('AGEB_aggregated.mat');

% Table census includes census data organized by AGEB ID. We now load a
% file containing the location (centroid) and Area of all AGEBS in GDL
% (coordinates are EPSG 32613 and areas are square meters)

geography = readtable('ageb_area_xy.csv');

geography = geography(:,[2,8,9,10]); %AGEB ID, AREA, X, Y

%We now seek to expand this data with demography from the census table

ID.geo = table2array(geography(:,1));
ID.cen = table2array(census(:,1));
ID.index = nan(size(geography,1),1);

for k = 1:size(geography,1)
    
    I = find(ID.cen == ID.geo(k));
    ID.index(k) = I;
    
end

geography = [geography, census(ID.index,2:end)];

demog = NaN(size(geography,1),2);

%Now replace missing data (marked as '*') for NaNs
for k = 1:size(geography,1)
    
    if string(geography.psinder(k)) == '*'
        demog(k,2) = nan;
    else
        demog(k,2) = str2num(cell2mat(geography.psinder(k)));
    end
    if string(geography.p_15ymas(k)) == '*'
        demog(k,1) = nan;
    else
        demog(k,1) = str2num(cell2mat(geography.p_15ymas(k)));
    end
    
end

%We now have a table with AGEB ID, Area, X,Y, Total population, population
%15 or older, and uninsured population

census = [geography(:,1:5), array2table(demog,'VariableNames',{'p_15ymas','psinder'})];

clear demog geography k I ID 

save('AGEB_aggregated.mat');


%% INTERPOLATE CENSUS DATA ONTO NUMERICAL GRID

load('AGEB_aggregated.mat'); %load census data and grid
load('/home/noel/Documents/RUDDS Final/RUDDS_Software/malla892nod.mat'); 

%Data structure grid includes variablas "mm" and "z". "mm" makes references to x
%and y components in "z" that make up each of all 1580 triangular elements
%in our numerical. 

%Start by calculating the centroid of all elements.
for k = 1:grid.nel
    
    grid.centroids(k,1) = mean(grid.z(1,grid.mm(:,k))); %x coordinate
    grid.centroids(k,2) = mean(grid.z(2,grid.mm(:,k))); %y coordinate
    
end

%Now interpolate data from census table onto these locations
city.xy = grid.centroids;


%Calculate the area of grid elements
for k = 1:grid.nel
    
    city.area(k) = 0.5*abs( grid.z(1,grid.mm(1,k))*(grid.z(2,grid.mm(2,k)) - ...
        grid.z(2,grid.mm(3,k))) + grid.z(1,grid.mm(2,k))*(grid.z(2, ...
        grid.mm(3,k)) - grid.z(2, grid.mm(1,k))) + grid.z(1, grid.mm(3,k))* ...
        (grid.z(2,grid.mm(1,k)) - grid.z(2,grid.mm(2,k))));
    
end

city.area = city.area .*1e-6;


%Exclude locations with 0 inhabitants or no data
pobcheck = census.pobtot ~= 0;

city.density = griddata(census.x(pobcheck), census.y(pobcheck), ...
    census.pobtot(pobcheck)./(census.area(pobcheck)/1e6),...
    city.xy(:,1), city.xy(:,2));

nancheck = ~isnan(census.p_15ymas);
nancheck = logical(nancheck.*pobcheck);

city.childratio = griddata(census.x(nancheck), census.y(nancheck), ...
    1-census.p_15ymas(nancheck)./census.pobtot(nancheck), city.xy(:,1), ...
    city.xy(:,2)); %fraction of people below 15 years old

nancheck = ~isnan(census.psinder);
nancheck = logical(nancheck.*pobcheck);

city.healthratio = griddata(census.x(nancheck), census.y(nancheck), ...
    census.psinder(nancheck)./census.pobtot(nancheck),city.xy(:,1), ...
    city.xy(:,2)); %fraction of people who don't have health insurance

%Using the density variable, we'll integrate over grid element area
%to get the total population at each location

city.population = city.area' .* city.density;

%%Visualize a variable
%load('/home/noel/Documents/MATLAB/Functions/my_colormaps.mat');

% 
% figure; hold on;
% fill(X,Y,city.childratio);
% 
% colormap(my_colormaps.goodjet);
% xlabel('Longitude (EPSG 32613)'); ylabel('Latitude (EPSG 32613)');
% xticks([6.6 6.7 6.8]*1e5); yticks([22.7 22.8 22.9 23]*1e5);



%Save city characteristics
clear grid census k nancheck pobcheck 

save('city.mat');

%Now run the script get_TAS to add some extra variables to our data
%structure city.mat



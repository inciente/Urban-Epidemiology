clc
clear all


%In this script, we take in data from the 2017 National Statistic Directory
%of Economical Units (DENUE) and the National System for School Information
%(SNIE) to calculate adult and child TAS at all our grid elements.

%Before running this, make sure you executed onto_grid.m and produced the
%file city.mat with no issues.

%Send your questions to:
% Noel Brizuela | nogutier@ucsd.edu
% Scripps Institution of Oceanography, University of California, San Diego
% December of 2018


%% DENUE

DENUE = readtable('denue_inegi_2016_redux.csv');
DENUE.latitud = str2double(table2array(DENUE(:,3)));
DENUE.codigo_act = str2double(table2array(DENUE(:,1)));
DENUE.act_redux = floor(table2array(DENUE(:,1))/10000);

%Now, transform the data coordinates onto EPSG 32613
DENUE.longitud = DENUE.longitud*1.0597e5 + 1.16238e7;
DENUE.latitud = DENUE.latitud*1.1007e5 + 1.11942e4;

%These are state-wide data, so we define a square and erase data outside 
X = [6.55e5, 6.55e5, 6.9e5, 6.9e5];
Y = [2.26e6, 2.3e6, 2.3e6, 2.26e6];

in = inpolygon(DENUE.longitud, DENUE.latitud, X,Y);
DENUE = DENUE(in,:);

%2 things left to do: convert columns 1 (Type of establishment) and 2
%(Employee number range) into TAS-relevant quantities. 

%Begin with categorizing businesses into Primary, Industry, Retail and
%Service. We load TAS_activity_codes.csv to establish equivalences. Modify 
%this file if, for example, you want to evaluate Hotels and restaurants as
%retail instead of services.

TAS_code = readtable('TAS_activity_codes.csv');

%Assign TAS coefficients for [Retail, Services, Industry, Primary, Student]
TAS_coef = [3.0, 2.4, 1.9, 1.3, 1.3];
TAS_cat = ["Retail","Services","Industry","Primary","Student"]; %students must be at the end

TAS_code = [TAS_code, table(NaN(size(TAS_code,1),1))];
%TAS_code.TAS_rate = table(zeros(size(TAS_code,1),1));
TAS_code.Properties.VariableNames{'Var1'} = 'TAS_rate';

for k = 1:size(TAS_code,1) 
    check = string(table2array(TAS_code(k,3))) == TAS_cat;
    TAS_code(k,4) = table(TAS_coef(check));
end

%Now translate text to number of workers for each entry in DENUE,
%cross-reference to TAS_code and calculate number of trips induced by each
%establishment

%What kind of establishments are these?
check = table2array(DENUE(:,5)) ==  table2array(TAS_code(:,1))';

%based on answer, asign TAS_rate to each one
DENUE.rate = sum(check.*TAS_code.TAS_rate',2);

%Now assign number of workers and multiply times rate
work_num = [3, 8, 20, 40, 75, 175, 290];
work_cat = ["0 a 5 personas","6 a 10 personas","11 a 30 personas", ...
    "31 a 50 personas","51 a 100 personas","101 a 250 personas", ...
    "251 y más personas"];

check = string(table2array(DENUE(:,2))) == work_cat;
DENUE.trips = sum(check.*work_num,2); %number of workers

%Multiply times rate to get number of trips
DENUE.trips = DENUE.trips .* DENUE.rate;

clear check work_num work_cat TAS_cat

%Now we assign each business to a grid element and get the total TAS from
%job sources. Still need to add adult students and calculate infantile TAS

load('malla892nod.mat');
load('city.mat')

%Go element by element 
for k = 1:length(grid.X)
   
    in = inpolygon(DENUE.longitud, DENUE.latitud, grid.X(:,k), grid.Y(:,k));
    city.TAS_adult(k) = sum(DENUE.trips(in));
    
end

city.TAS_adult = city.TAS_adult'; clear in k 

%Great! we're almost done with DENUE. Now just save the locations of
%private and state-owned hospitals (SCIAN CODE = 622111, 622112)

hosp = table2array(DENUE(:,1)) == 622111;
city.hospitals.private = [table2array(DENUE(hosp,4)), ...
    table2array(DENUE(hosp,3))];

hosp = table2array(DENUE(:,1)) == 622112;
city.hospitals.public = [table2array(DENUE(hosp,4)), table2array(DENUE(hosp,3))];

clear hosp TAS_code


%% SNIE

high_school = load('prepasGDL.csv');
university = load('unisGDL.csv');

%Data array is # of students, longitude, latitude
lower_ed = load('escuelasGDL.csv');
higher_ed = cat(1,high_school, university);

clear high_school university

%Transform coordinates to EPSG 32613. Higher ed is offset by ~ 1.2 km west
lower_ed(:,2) = lower_ed(:,2)*1.0597e5 + 1.16238e7;
lower_ed(:,3) = lower_ed(:,3)*1.1007e5 + 1.11942e4;

higher_ed(:,2) = higher_ed(:,2)*1.05967e5 + 1.16250e7;
higher_ed(:,3) = higher_ed(:,3)*1.1007e5 + 1.11942e4;

%These are state-wide data, so we define a square and erase data outside 
X = [6.55e5, 6.55e5, 6.9e5, 6.9e5];
Y = [2.26e6, 2.3e6, 2.3e6, 2.26e6];

in = inpolygon(lower_ed(:,2), lower_ed(:,3), X,Y);
lower_ed = lower_ed(in,:);

in = inpolygon(higher_ed(:,2), higher_ed(:,3), X,Y);
higher_ed = higher_ed(in,:);

%Now go through all grid elements and add their TAS

for k = 1:length(grid.X)
    
    in = inpolygon(lower_ed(:,2), lower_ed(:,3), grid.X(:,k), grid.Y(:,k));
    city.TAS_child(k) = sum( lower_ed(in,1))*TAS_coef(end);
    
    in = inpolygon(higher_ed(:,2), higher_ed(:,3), grid.X(:,k), grid.Y(:,k));
    city.TAS_adult(k) = city.TAS_adult(k) + sum(higher_ed(in,1))*TAS_coef(end);
    
end

%Data structure city that results after this step is what we use as input
%to calculate origin-destination matrices and run our epidemiological
%simulations.

city.TAS_child = city.TAS_child';
city.area = city.area';

city.child = city.population .* city.childratio;
city.adult = city.population - city.child;

save('city.mat','city');

%% Visualize TAS data sources

% load('malla892nod.mat');

figure; ax(1) = subplot(131); hold on;
fill(grid.X,grid.Y,NaN(1580,1));
scatter(DENUE.longitud, DENUE.latitud,'.','m');
title('Jobs')


ax(2) =subplot(132); hold on;
fill(grid.X, grid.Y, NaN(1580,1));
scatter(higher_ed(:,2), higher_ed(:,3),'b','.');
title('Higher education')

ax(3) = subplot(133); hold on;
fill(grid.X, grid.Y, NaN(1580,1));
scatter(lower_ed(:,2), lower_ed(:,3), 'k','.');
title('Lower education')

linkaxes(ax);

for k = 1:3
    subplot(1,3,k)
    xticks([6.6 6.7 6.8]*1e5); 
    if k == 1
        yticks([22.7 22.8 22.9 23]*1e5);
    end
end

a= -103.3272*1.0597e5 + 1.16238e7;
b= 20.6864*1.1007e5 +1.11942e4;  











clc
clear all

% Master script to setup and run epidemiological simulations. 
% You will need to load a city.mat file (you can generate your own using
% the scripts in the Demography folder) or create an analogue with the
% variables for your region of interest.


%% SET UP DISEASE PARAMETERS

infected.beta.adult = 1.25e-5;
infected.beta.child = 3.47e-5;
infected.gamma = 1/4;

infected.isolation = 0.8; %coefficient alpha to reduce infected mobility
infected.vishosp = 0.08; %probability of visit to hospital during infective period
infected.bed = 0.03; %probability of being hospitalized after visiting (for calibration purposes)

infected.t0 = 15; %number of infecteds at time of onset
infected.age = 0; %is patient zero an adult (0) or a child (1)?

in_cond = 1; %number of initial conditions you want to test 
ndays = 250; %number of days you wish to simulate

addpath('/home/noel/Documents/MATLAB/Functions/')
addpath('/home/noel/Documents/RUDDS Final/RUDDS_Software/Demography')
load('city.mat');
load('malla892nod.mat');
load('mobility.mat');


offset = 45; % time offset (days) to plot between simulation results and data

%% GET MOBILITY - COMMENT OUT IF YOU ALREADY HAVE OD MATRICES

clear mob

%First of all, we create a matrix of distances between centroids of all our
%grid elements.
r_xy = sqrt( (city.xy(:,1) - city.xy(:,1)').^2 + ...
    (city.xy(:,2) - city.xy(:,2)').^2)/1000;


%Now, select the mobility parameters that define P(\Delta r_xy)
mob.r_0 = 5; 
mob.b = 1.75;
mob.k = 80;

%Call get_OD_matrix to calculate the origin-destination matrices of
%susceptibles.
mob.ODadult = get_OD_matrix(city.TAS_adult',r_xy,mob.r_0,mob.b,mob.k);
mob.ODchild = get_OD_matrix(city.TAS_child',r_xy,mob.r_0,mob.b,mob.k);

%Multiply by sum(TAS)/sum(population) to adjust for number of daily trips
mob.ODchild = mob.ODchild * (sum(city.TAS_child)/sum(city.child));
mob.ODadult = mob.ODadult * (sum(city.TAS_adult)/sum(city.adult));

%locate each hospital in its corresponding grid element
for k = 1:length(grid.X)
    
    in = inpolygon(city.hospitals.public(:,1), city.hospitals.public(:,2), ...
        grid.X(:,k), grid.Y(:,k));   
    if sum(in) == 1
        city.hospitals.public(in,3) = k;
    end
    
end

%Delete hospitals outside our numerical grid.
city.hospitals.public(city.hospitals.public(:,3) == 0,3) = NaN;
nancheck = ~isnan(city.hospitals.public(:,3));
city.hospitals.correspondence = city.hospitals.public(nancheck,3);

%Calculate distance from grid elements to hospitals.
r_xy = sqrt( (city.xy(:,1) - city.xy(city.hospitals.correspondence,1)').^2 + ...
    (city.xy(:,2) - city.xy(city.hospitals.correspondence,2)').^2 );

[dist, corresp] = min(r_xy'); %corresp tells us which hospital people from 
%each grid element will attend. 

corresp = city.hospitals.correspondence(corresp);
mob.ODhosp = zeros(size(mob.ODadult));

clear city.hospitals.correspondence;

%For each origin, add 1 to the destination of nearest hospital
for k = 1:size(mob.ODhosp,1)  
    hospital = corresp(k);
    mob.ODhosp(k,hospital) = 1;
end

clear r_xy nancheck k hospital  dist corresp in

save('/home/noel/Documents/RUDDS Final/RUDDS_Software/mobility.mat','mob');


%% RUN SIMULATION 

%Adjust value of beta to variability in area of grid elements
infected.beta.adult = infected.beta.adult./city.area;
infected.beta.child = infected.beta.child./city.area;

%Get the OD matrices of infected people
infected.ODadult = infected.isolation*mob.ODadult + ...
    infected.vishosp*infected.gamma*mob.ODhosp; 

infected.ODchild = infected.isolation*mob.ODchild + ...
    infected.vishosp*infected.gamma*mob.ODhosp;

%Save OD matrices in GPU for faster computing 
infected.ODadult = gpuArray(infected.ODadult);
infected.ODchild = gpuArray(infected.ODchild);

mob.ODadult = gpuArray(mob.ODadult);
mob.ODchild = gpuArray(mob.ODchild);

% Compute R_0

%Let's first see how susceptibles are distributed across the city:
susc_city = (mob.ODadult')*city.adult + (mob.ODchild')*city.child;

%Now multiply by probability of sick person visiting those places:
R_0.adult = (infected.ODadult*(susc_city.*infected.beta.adult))./infected.gamma; 
R_0.child = (infected.ODchild*(susc_city.*infected.beta.child))/infected.gamma;

%Get weighted average of R_0 using population statistics 
R_0.mean_adult = gather(sum(R_0.adult.*city.adult)/sum(city.adult));
R_0.mean_child = gather(sum(R_0.child.*city.child)/sum(city.child));


% Create arrays to store solution
Solution = NaN(in_cond, ndays); 
c_cases = NaN(in_cond, ndays+1);
c_cases(:,1) = infected.t0; 



%New Hospitalizations: %First column for adults, second for children.
hosp_new = NaN(in_cond, ndays,2); 
attack_rate_adult = zeros(1580, in_cond); 
attack_rate_child = zeros(1580, in_cond);

tic
%Begin simulation and cycle through all possible origins of epidemic
for origin=1:in_cond
    
    
    %HOSPITALIZATION BUFFER: Hospitalizations will occur until the third
    %day of illness, so that infectives from day 1 are counted as hospitalized
    %until day 3: First column is adults, second column is children.
    hosp_buffer = zeros(3,2);
    hosp_new2 = NaN(ndays,2);
    
    
    %Create arrays to store number of S, I, R at all parts of the city.
    SIR_Adult = zeros(size(city.xy,1),3); 
    SIR_Adult = gpuArray(SIR_Adult); 
    
    SIR_Children = zeros(size(city.xy,1),3); 
    SIR_Children = gpuArray(SIR_Children);

    %Condition when everyone is susceptible:
    SIR_Adult(:,1) = city.adult; 
    SIR_Children(:,1) = city.child; 
    
    %Set the starting condition for I(x,t=0) 

    if infected.age == 1
        SIR_Children(origin,2) = infected.t0;
        SIR_Children(origin,1) = SIR_Children(origin,1) - infected.t0;
    else
        SIR_Adult(origin,2) = infected.t0; 
        SIR_Adult(origin,1) = SIR_Adult(origin,1) - infected.t0; 
    end
    
    %Now the actual solution for a given origin:
        
    for time = 1:ndays
        
        %Compute how many infected people visit each point of the city
        I_map = (infected.ODadult')*SIR_Adult(:,2).*infected.beta.adult  + ...
            (infected.ODchild')*SIR_Children(:,2).*infected.beta.child; 
        
        %Compute expected number of susceptibles visiting each place (need
        %to store where they come from, so we just do an element-wise
        %multiplication of OD matrix and array S.
        
        S_trip_Adult = mob.ODadult.*SIR_Adult(:,1);
        S_trip_Children = mob.ODchild.*SIR_Children(:,1);
        
        %Change in number of infected individuals over this timestep
        di_Adult = S_trip_Adult*(I_map); 
        di_Children = S_trip_Children*(I_map); 

        %Number of recovered people: 
        SIR_Adult(:,3) = SIR_Adult(:,3) + infected.gamma*SIR_Adult(:,2); 
        SIR_Children(:,3) = SIR_Children(:,3) + infected.gamma*SIR_Children(:,2);
        
        SIR_Adult(:,1) = SIR_Adult(:,1) - di_Adult; 
        SIR_Adult(:,2) = SIR_Adult(:,2) + di_Adult - ...
            infected.gamma*SIR_Adult(:,2); 
        
        SIR_Children(:,1) = SIR_Children(:,1) - di_Children; 
        SIR_Children(:,2) = SIR_Children(:,2) + di_Children - ...
            infected.gamma*SIR_Children(:,2); 
       
        
        SIR_Children(SIR_Children < 0) = 0; 
        SIR_Adult(SIR_Adult < 0) = 0;
    

        %Save I(x,t) to solution file. -
        Solution(origin, time) = gather(sum(SIR_Adult(:,2) + SIR_Children(:,2)));
        
        %CUmulative cases since beginning of epidemic
         c_cases(origin, time+1) = gather(sum(di_Adult) + sum(di_Children) ...
             + c_cases(origin, time));
        
       attack_rate_adult(:,origin) = attack_rate_adult(:,origin) + ...
           gather(di_Adult)./city.adult;
       
       attack_rate_child(:,origin) = attack_rate_child(:,origin) + ...
           gather(di_Children)./city.child;
        
        
        
        %Save number of new infectives to compute new hospitalizations
        hosp_buffer(1,1) = gather(sum(di_Adult)); 
        hosp_buffer(1,2) = gather(sum(di_Children)); 
        
        %Use new infectives from time-2 to get new hospitalizations.
        hosp_new2(time, 1) = gather(hosp_buffer(3,1)* ...
            infected.vishosp*infected.bed);
        hosp_new2(time, 2) = gather(hosp_buffer(3,2)* ...
            infected.vishosp*infected.bed);
        
        %Make room for numbers from next timestep
        hosp_buffer(2,:) = hosp_buffer(1,:); 
        hosp_buffer(3,:) = hosp_buffer(2,:); 
        

        
    end
    
    hosp_new(origin,:,:) = hosp_new2;
    
end

toc

clear susc_city I_map S_trip_Adult S_trip_Children 
clear di_Adult di_Children 


%% COMPARE RESULTS TO 2009 DATA



load AH1N1_data.mat

figure(1);

set(0,'DefaultAxesFontsize',18);


hold on;

for i=1:in_cond
    
    
curveHosp=hosp_new(i,:,1)+hosp_new(i,:,2);

plot(linspace(1,length(curveHosp),length(curveHosp)),curveHosp, ...
    'Color',[1 1 1]*0.45)
hold on

end

line1=plot((7:7:33*7) - offset,data,'bx');
set(line1,'markerSize',12,'LineWidth',2);

axis([0 33*7+20 0 60])

xlabel('Time (days)')
ylabel('New hospitalizations')

%save('results-calibration_Feb18.mat','-mat')



R_0



















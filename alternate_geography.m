function epidemic = alternate_geography(spatial_TAS, spatial_DENS, expop, extrip, in_cond)


% Run large amounts of epidemiological simulations varying the spatial
% structure and inequalities of a city. 
% spatial_TAS = 0   to use spatial distribution of TAS from city.mat
%             = 1   to randomnize spatial distribution


% spatial_DENS= 0   to use spatial distribution of people from city.mat
%             = 1   to randomnize distribution of people
%             = 2   to sort density as TAS


% expop and extrip are used to re-allocate housing and trip destinations
% respectively. values < 1 reduce Gini coefficient, 1 keeps it constant,
% and values > 1 increase Gini coefficient.

%in_cond is the number of initial conditions you want to test for each
%configuration. I suggest somewhere between 10 and 50, unless you're not
%concerned about running time. This only impacts our estimate of the
%attack rate, which tends to be stable anyway.

%% SET UP DISEASE PARAMETERS

infected.beta.adult = 1.25e-5;
infected.beta.child = 3.47e-5;
infected.gamma = 1/4;

infected.isolation = 0.8; %coefficient alpha to reduce infected mobility
infected.vishosp = 0.08; %probability of visit to hospital during infective period
infected.bed = 0.02; %probability of being hospitalized after visiting (for calibration purposes)

infected.t0 = 15; %number of infecteds at time of onset
infected.age = 1; %is patient zero an adult (0) or a child (1)?

ndays = 350; %number of days you wish to simulate


load('city.mat');
load('malla892nod.mat');


%% MODIFY GEOGRAPHY

pobtot.child  = sum(city.child);
pobtot.adult  = sum(city.adult);

%Modify the spatial structure of the city?

if spatial_TAS == 1
    
    %Randomnize TAS
    city.TAS_adult = city.TAS_adult(randperm(length(city.xy)));
    city.TAS_child = city.TAS_child(randperm(length(city.xy)));   
    
end

if spatial_DENS == 1
    
    %Randomnize density and distribution of children.
    city.density = city.density(randperm(length(city.xy)));
    city.childratio = city.childratio(randperm(length(city.xy)));
 
    city.adult = (1 - city.childratio).*city.density.*city.area;
    city.child = city.childratio.*city.density.*city.area;
    
    %Preserve size of populations
    city.adult = city.adult*pobtot.adult/sum(city.adult);
    city.child = city.child*pobtot.child/sum(city.child);
    
    %Update other variables 
    city.population = city.adult + city.child;
    city.density = city.population./ city.area; 
    
end

    
if spatial_DENS == 2
    
    [dummy, sort_tas_adult] = sort(city.TAS_adult);
    [dummy, sort_tas_child] = sort(city.TAS_child);
    
    [dummy, sort_dens_adult] = sort(city.adult./city.area);
    [dummy, sort_dens_child] = sort(city.child./city.area);
    
    %Give high adult density to neighborhoods with high adult TAS
    city.adult(sort_tas_adult) = city.adult(sort_dens_adult)./city.area(sort_dens_adult);
    city.adult = city.adult.*city.area;
    
    city.child(sort_tas_child) = city.child(sort_dens_child)./city.area(sort_dens_child);
    city.child = city.child.*city.area;
    
    %Preserve size of population
    city.adult = city.adult * pobtot.adult/sum(city.adult);
    city.child = city.child * pobtot.child/sum(city.child);
    
end

if extrip ~= 1
    % Change the degrees of inequality:
    city.TAS_adult = rescale_demo(city.TAS_adult,extrip);
    city.TAS_child = rescale_demo(city.TAS_child, extrip);
end

if expop ~= 1
    city.child = rescale_demo(city.child, expop);
    city.adult = rescale_demo(city.adult, expop);
end


city.population = city.child + city.adult;
city.density = city.population./city.area;
city.childratio = city.child./city.population;


%% GET MOBILITY - COMMENT OUT IF YOU ALREADY HAVE OD MATRICES

%%% clear mob

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


%% RUN SIMULATION 

%Adjust value of beta to variability in area of grid elements
infected.beta.adult = infected.beta.adult./city.area';
infected.beta.child = infected.beta.child./city.area';

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

% %Calculate second-generation R_0
% susc_city = (mob.ODadult')*(city.adult.*R_0.adult) + ...
%     (mob.ODchild')*(city.child.*R_0.child);
% 
% R_0.dob_adult = (infected.ODadult*(susc_city.*infected.beta.adult))./infected.gamma; 
% R_0.dob_child = (infected.ODchild*(susc_city.*infected.beta.child))/infected.gamma;


% Create arrays to store solution
Solution = NaN(in_cond, ndays); 
c_cases = NaN(in_cond, ndays+1);
c_cases(:,1) = infected.t0; 



%New Hospitalizations: %First column for adults, second for children.
hosp_new = NaN(in_cond, ndays,2); 
%attack_rate_adult = zeros(1580, in_cond); 
%attack_rate_child = zeros(1580, in_cond);


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
        
        SIR_Adult(:,1) = SIR_Adult(:,1) - di_Adult; 
        SIR_Adult(:,2) = SIR_Adult(:,2) + di_Adult - ...
            infected.gamma*SIR_Adult(:,2); 
        
        SIR_Children(:,1) = SIR_Children(:,1) - di_Children; 
        SIR_Children(:,2) = SIR_Children(:,2) + di_Children - ...
            infected.gamma*SIR_Children(:,2); 
       
        %Number of recovered people: 
        SIR_Adult(:,3) = SIR_Adult(:,3) + infected.gamma*SIR_Adult(:,2); 
        SIR_Children(:,3) = SIR_Children(:,3) + infected.gamma*SIR_Children(:,2);
        
        SIR_Children(SIR_Children < 0) = 0; 
        SIR_Adult(SIR_Adult < 0) = 0;
    

        %Save I(x,t) to solution file. -
        Solution(origin, time) = gather(sum(SIR_Adult(:,2) + SIR_Children(:,2)));
        
        %CUmulative cases since beginning of epidemic
        c_cases(origin, time+1) = gather(sum(di_Adult) + sum(di_Children) ...
            + c_cases(origin, time));
        
%        attack_rate_adult(:,origin) = attack_rate_adult(:,origin) + ...
%            gather(di_Adult)./city.adult;
%        
%        attack_rate_child(:,origin) = attack_rate_child(:,origin) + ...
%            gather(di_Children)./city.child;
        
        
        
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
 


clear susc_city I_map S_trip_Adult S_trip_Children 
clear di_Adult di_Children 


gini_o = Gini_curve(city.population, city.area);

work_hours = mob.ODadult'*city.adult + mob.ODchild'*city.child;
gini_d = Gini_curve(work_hours, city.area);

epidemic.R_0_child = gather([min(R_0.child), R_0.mean_child, max(R_0.child)]);
epidemic.R_0_adult = gather([min(R_0.adult), R_0.mean_adult, max(R_0.adult)]);

epidemic.attack_rate = mean(c_cases(:,end))./sum(city.population);

epidemic.gini_o = gini_o.coeff;
epidemic.gini_d = gather(gini_d.coeff);



















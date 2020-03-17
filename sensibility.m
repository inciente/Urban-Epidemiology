clc
clear all

this_path = pwd;
addpath([this_path,'/Demography']);

% Script to run lots of epidemiological simulations in alternate
% geographies. Read the script alternate_geography.m to learn more about
% how this is achieved. 
repetition = 1;

% LIST OF EXPONENTS TO TRY
extrip = linspace(0.4, 1.3,2);
expop = linspace(0.7, 1.3, 2);

extrip = repmat(extrip,1, repetition); %repeat many times for averaging?
expop = repmat(expop, 1, 1); 


gini_o = NaN(length(extrip), length(expop));
gini_d = NaN(length(extrip), length(expop));

attack_rate = NaN(length(extrip), length(expop));

%store minimum, mean and maximum (3) 
R_0.child = NaN(length(extrip), length(expop),3);
R_0.adult = NaN(length(extrip), length(expop),3);



for k = 1:length(extrip)
    tic
    for m = 1:length(expop)
    
        epid = alternate_geography(1, 1, expop(m), extrip(k), 5);
        
        gini_o(k,m) = epid.gini_o;
        gini_d(k,m) = epid.gini_d;
        
        attack_rate(k,m) = epid.attack_rate;
        
        R_0.child(k,m,:) = epid.R_0_child;
        R_0.adult(k,m,:) = epid.R_0_adult;
        
        
    end
    toc
        
end



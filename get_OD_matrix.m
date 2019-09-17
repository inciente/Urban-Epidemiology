function OD_matrix = get_OD_matrix(TAS, r_xy, r_0, b, k)

% Function to calculate probabilistic origin-destination matrices using the
% methodology described in 'Understanding the Role of Urban Design in
% Disease Spreading'. 
%
% Input: 
%
% TAS     -  Trip Attraction Strength for grid elements in a city grid. See 
%            scripts onto_grid.m and get_TAS.m for more details. 
% r_xy    -  Distance matrix that gives the distance (in km) between the
%            centroids of grid elements. 
% r_0,b,k -  are all parameters that define the shape of the 'willingness
%            to travel' curve. These can be scalars OR vectors of the same
%            dimensions as TAS

% Output:
% OD_matrix  This is a matrix of the same dimension as r_xy that gives us
%            the probability that a susceptible visits each of the
%            different grid elements on a given day. It is a function of
%            mobility parameters, so a different one of these will have to
%            calculated for every age group taken into account. This matrix
%            will be modified elsewhere to infer the mobility of infecteds.
%            Dimensions are (origin, destination)
%
%            We assume that c(x) = 1 for all locations, but this can be
%            changed by multiplying the whole matrix or row to different
%            values of c(x). 
%
%
% Send your questions to:
% Noel Brizuela | nogutier@ucsd.edu
% Scripps Institution of Oceanography. University of California, San Diego
% December of 2018


if isequal(size(r_0),size(b)) && isequal(size(b),size(k))
    % proceed
else
    warning('All mobility parameters r_0, b, k must have the same dimension and orientation')
    return
end

if size(r_0,1) < size(r_0,2)
    r_0 = r_0';
    b = b';
    k = k';
end

if size(TAS,1) > size(TAS,2)
    TAS = TAS';
end

OD_matrix = (r_xy + r_0).^(-b).*exp(-r_xy./k).*TAS;

%Normalize probability to c = 1
OD_matrix = OD_matrix./sum(OD_matrix,2);

end
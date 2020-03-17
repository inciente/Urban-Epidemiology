function [new_var] = rescale_demo(variable, exponent)

%Function to rescale demographic variables using a given exponent

new_var = (variable/(mean(variable))).^exponent; 
new_var = new_var / mean(new_var) * mean(variable); 

end

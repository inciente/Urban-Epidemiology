function lorenz = Gini(variable, area)

%Function to obtain the lorenz curve and Gini coefficient for a geographic 
%variable.


if size(variable,1) == size(area,2)
    area = area';
end


[new_var, index] = sort(variable./area);

lorenz.y = cumsum(new_var)/sum(new_var); 
lorenz.x = cumsum(area(index))/sum(area);

percent = cumsum(ones(length(variable),1))/length(variable); 

%Area below data
lorenz.coeff = sum(lorenz.y.*area(index))/sum(area);

%Area below line of equality is 0.5 so:
lorenz.coeff = (0.5 - lorenz.coeff)/0.5;

end

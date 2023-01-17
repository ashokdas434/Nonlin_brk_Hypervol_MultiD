%% This function creates the grid and pivot points
function [x,R,del_x] = Grids(x_min,x_max, I)
%%
x = zeros(1,I); del_x = zeros(1,I); % Initialization

if x_min ==0
    x_min = 1e-6;
end
%%
R =  exp(linspace(log(x_min),log(x_max),I+1));

if x_min ==0
    R(1) = 0;
end
%%
for i=1:I
    x(i)     = 0.5*(R(i)+R(i+1));
    del_x(i) = R(i+1) - R(i);
end

return
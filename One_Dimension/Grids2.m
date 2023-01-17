%% This function creates the grid and pivot points (Including 1 as pivot)
function [x,R,del_x, L] = Grids2(x_min,x_max, I, grid_mul)
%%
x = zeros(1,I); del_x = zeros(1,I); % Initialization

if x_min ==0
    x_min = 1e-8;
end
%%
R =  exp(linspace(log(x_min),log(x_max),I+1));

if x_min ==0
    R(1) = 0;
end

%% Division into smaller cells using multiplier
if grid_mul >1
    R1 = R(1);
    for i=1:I
        del = R(i+1)-R(i);
        for j=1:grid_mul
            R1 = [R1, R1(end)+del/grid_mul];
        end
    end
    R = R1;
    L = length(R)-1;  % No. of cells
else
    L = I;
end
%%
R(L+1) = R(L+1) + (R(L+1)-R(L));
%%
for i=1:L
    x(i)     = 0.5*(R(i)+R(i+1));
    del_x(i) = R(i+1) - R(i);
end

return
%% This function creates the matrix version of the collision function
function K = K_Fun(example,x,I)
K = zeros(I); % initialization

switch example
    case 2 % K = 1
        K = ones(I);
        
    case 1 % K_st(i,j) = i*j
        for i=1:I
            j=1:I;
            K(i,j)= x(i)* x(j);
        end
end

return
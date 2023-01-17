%% This function creates the matrix version of the collision function
function K = K_Fun(K_index,x1,x2,I1,I2)
K = zeros(I1,I2,I1,I2); % initialization

switch K_index
    case 1 % K = 1
        K = ones(I1,I2,I1,I2);
        
    case 2 % K_st(i,j) = i+j
        for i=1:I1
            for j=1:I2
                k=(1:I1);
                l=1:I2;
                K(i,j,:,:)= (x1(i)*x2(j))*(x1(k)'*x2(l));
            end
        end
end

return
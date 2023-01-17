%% Function to formulate the weights
function [w1,w2_b,w2_d] = weights(x1,x2,B)

I1 = length(x1); I2 = length(x2);
w1 = zeros(I1,I2); w2_b = zeros(I1,I2); w2_d = zeros(I1,I2); % Initialization
nu = 4; % number of fragments per breakage

x1x2 = x1'*x2;

for i=1:I1
    for j=1:I2
        m=1:i; n=1:j;
        S = x1(m) * (B(m,n,i,j)*x2(n)');
        w1(i,j) = S/(x1(i)*x2(j));
        
        temp = x1(i)*x2(j)- x1x2(m,n);
        temp2 = temp.* B(m,n,i,j);
        S2 = sum(sum(temp2));
        w2_b(i,j) = x1(i)*x2(j)* (nu -1)/S2;
        w2_d(i,j) = w2_b(i,j) * w1(i,j);
    end
end

w2_b(1,1) = 0; w2_d(1,1) = 0;

return
    

%% Function to calculate F for consevative approach
function F = F_conserve(x,del_x,N,K,beta)

I = length(x); F = zeros(1,I+1); % Initialization

x_del_x = x.*del_x;
KN  = K*N; % Multiplication of K & N (I*1 vector)
NKN = N.*KN;

for i=1:I
    x_del_x_beta = x_del_x(1:i) * beta(1:i,i+1:I);
    
    F(i+1) = x_del_x_beta * NKN(i+1:I);
end

return

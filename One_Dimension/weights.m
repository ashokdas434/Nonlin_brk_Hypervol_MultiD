%% Function to formulate the weights
function [w1,w2_b,w2_d] = weights(x,B,frag)

I = length(x);
w1 = zeros(1,I); w2_b = zeros(1,I); w2_d = zeros(1,I); % Initialization
%%
for i=1:I
    j = 1:i;
    w1(i) = ( x(j)*B(j,i) )/x(i);
    
    w2_b(i) = x(i)*(frag-1)/((x(i)-x(j))*B(j,i));
    w2_d(i) = w2_b(i)*w1(i);
end

w2_b(1) = 0; w2_d(1) = 0;


return
    
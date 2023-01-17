%% Discrete rate function for Mass conserving (MC) technique
function dNdt = discrete_MC(t,N, K,B,w1,x) 

I = length(N);
dNdt = zeros(I,1);

KN  = K*N; % Multiplication of K & N (I*1 vector)
NKN = N.*KN;

for i=1:I
    j = i:I;
    dNdt(i) = B(i,j)*NKN(j) - N(i)*w1(i)*KN(i);
end
        


fprintf('FVS-MC| t_sim=%1.1f | t_real=%1.8f | N_p=%1.4f | M=%1.5f\n',...
    toc, t, sum(N), x*N)

return
%% Discrete rate function for Number preserving & Mass conserving (NPMC) technique
function dNdt = discrete_NPMC(t,N, K,B,w2_b,w2_d,x) 

I = length(N);
dNdt = zeros(I,1);

KN  = K*N; % Multiplication of K & N (I*1 vector)
NKN = N.*KN;

for i=1:I
    j = i:I;
    wbB = w2_b.* B(i,:);
    dNdt(i) = wbB(j)*NKN(j) - N(i)*w2_d(i)*KN(i);
end

fprintf('FVS-NPMC| t_sim=%1.1f | t_real=%1.8f | N_p=%1.4f | M=%1.5f\n',...
    toc, t, sum(N), x*N)

return
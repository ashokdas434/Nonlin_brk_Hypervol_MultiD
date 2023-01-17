%% Discrete rate function for conservative approach
function dNdt = discrete_conserve(t,N, x,del_x,K,beta) 

I = length(N);
dNdt = zeros(I,1);

F = F_conserve(x,del_x,N,K,beta); % Function for conservative approach

for i=1:I
    dNdt(i) = (F(i+1)-F(i))/x(i);
end

fprintf('FVS-conserve| t_sim=%1.1f | t_real=%1.8f | N_p=%1.4f | M=%1.5f\n',...
    toc, t, sum(N), x*N)

return
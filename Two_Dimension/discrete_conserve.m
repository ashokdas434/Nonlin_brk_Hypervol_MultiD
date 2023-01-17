%% Discrete rate function for conservative approach
function dNdt = discrete_conserve(t,N, x1,x2,del_x1,del_x2,K,beta) 

I1 = length(x1); I2 = length(x2); 

% N in matrix form
N_mat = vec2mat(N,I1,I2);

%%
[F1,F2,F3] = F_conserve_check(x1,x2,del_x1,del_x2,N_mat,K,beta); % Functions for conservative approach
%%
for i=1:I1
    for j=1:I2
        dNdt_mat(i,j) = ( (F1(i+1,j)-F1(i,j))*del_x2(j) +(F2(i,j+1)-F2(i,j))*del_x1(i) -...
                      (F3(i+1,j+1)-F3(i+1,j)-F3(i,j+1)+F3(i,j)) )/(x1(i)*x2(j));
    end
end

%% matrix to vector form
dNdt = mat2vec(dNdt_mat);

fprintf('FVS-Conserve| t_sim=%1.1f | t_real=%1.9f | N_p=%1.4f | M=%1.5f\n',...
    toc, t, sum(sum(N_mat)), x1*(N_mat*x2'))


return
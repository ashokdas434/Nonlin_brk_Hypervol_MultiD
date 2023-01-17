%% Discrete rate function for Number preserving & Mass conserving (NPMC) technique
function dNdt = discrete_NPMC(t,N, K,B,w2_b,w2_d,x1,x2) 

I1 = length(x1); I2 = length(x2); 

% N in matrix form
N_mat = vec2mat(N,I1,I2);

%%
KN = zeros(I1,I2);
for p=1:I1
    for q=1:I2
        KN(p,q) = sum(sum(K(:,:,p,q).*N_mat));
    end
end
NKN = N_mat.*KN;

%%
for i=1:I1
    for j=1:I2
        birth = 0;
        for m=i:I1
            for n=j:I2
                birth = birth+ w2_b(m,n)*B(i,j,m,n)*NKN(m,n);
            end
        end
        dNdt_mat(i,j) = birth - w2_d(i,j)*NKN(i,j);
    end
end

%% matrix to vector form
dNdt = mat2vec(dNdt_mat);

fprintf('FVS-NPMC| t_sim=%1.1f | t_real=%1.9f | N_p=%1.4f | M=%1.5f\n',...
    toc, t, sum(sum(N_mat)), x1*(N_mat*x2'))

return
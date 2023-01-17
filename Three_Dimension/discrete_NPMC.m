%% Discrete rate function for Number preserving & Mass conserving (NPMC) technique
function dNdt = discrete_NPMC(t,N, K,B,w2_b,w2_d,x1,x2,x3) 

I1 = length(x1); I2 = length(x2); I3 = length(x3);
dNdt_mat = zeros(I1,I2,I3); 

% N in matrix form
N_mat = reshape(N,I1,I2,I3);

%%
KN = zeros(I1,I2,I3);
for p=1:I1
    for q=1:I2
        for r=1:I3
            KN(p,q,r) = sum(sum(sum(K(:,:,:,p,q,r).*N_mat)));
        end
    end
end
NKN = N_mat.*KN;

%%
Hyp_vol = 0; Hyp_vol_ch =0;

for i1=1:I1
    for i2=1:I2
        for i3=1:I3
            birth = 0;
            for m1=i1:I1
                for m2=i2:I2
                    for m3=i3:I3
                        birth = birth+ B(i1,i2,i3,m1,m2,m3)*NKN(m1,m2,m3)*w2_b(m1,m2,m3);
                    end
                end
            end
            dNdt_mat(i1,i2,i3) = birth - w2_d(i1,i2,i3)*NKN(i1,i2,i3);
            Hyp_vol = Hyp_vol + x1(i1)*x2(i2)*x3(i3)*N_mat(i1,i2,i3);
            Hyp_vol_ch = Hyp_vol_ch + x1(i1)*x2(i2)*x3(i3)*dNdt_mat(i1,i2,i3);
        end
    end
end

%% matrix to vector form
dNdt = reshape(dNdt_mat,I1*I2*I3,1);

fprintf('FVS-MC| t_sim=%1.1f | t_real=%1.5f | N_p=%1.4f | del_M= %1.5f | M=%1.5f\n',...
    toc, t, sum(sum(sum(N_mat))),Hyp_vol_ch, Hyp_vol)


return
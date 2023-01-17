% Function to calculate error functions
clear all
Example = 2; I = 30;
N_sim = 5;

for i = 1:N_sim
    G = 2^(i-1);
    load(['bEx-',num2str(Example),'-Grid_mul-',num2str(G),'.mat']);

    %    N_MC = N_MC./del_x; N_NPMC = N_NPMC./del_x; N_cons = N_cons./del_x; N_ANA = N_ANA./del_x;
    for j=1:I
        N_MC_all(i,j)   = sum(N_MC(G*(j-1)+1:G*j));
        N_NPMC_all(i,j) = sum(N_NPMC(G*(j-1)+1:G*j));
        N_cons_all(i,j) = sum(N_cons(G*(j-1)+1:G*j));
        N_ANA_all(i,j)  = sum(N_ANA(G*(j-1)+1:G*j));
    end
end

%% MC
for i=1:N_sim
    L1_MC(i) = sum(abs(N_MC_all(i,:)- N_ANA_all(i,:)))/sum(abs(N_ANA_all(i,:)));
    L2_MC(i) = sum((N_MC_all(i,:)- N_ANA_all(i,:)).^2)/sum((N_ANA_all(i,:)).^2);
    if i~=N_sim
        EOC_MC(i) = log(sum(abs(N_MC_all(i,:)- N_ANA_all(i,:)))/sum(abs(N_MC_all(i+1,:)- N_ANA_all(i+1,:))))/log(2);
    end
end

%% NPMC
for i=1:N_sim
    L1_NPMC(i) = sum(abs(N_NPMC_all(i,:)- N_ANA_all(i,:)))/sum(abs(N_ANA_all(i,:)));
    L2_NPMC(i) = sum((N_NPMC_all(i,:)- N_ANA_all(i,:)).^2)/sum((N_ANA_all(i,:)).^2);
    if i~=N_sim
        EOC_NPMC(i) = log(sum(abs(N_NPMC_all(i,:)- N_ANA_all(i,:)))/sum(abs(N_NPMC_all(i+1,:)- N_ANA_all(i+1,:))))/log(2);
    end
end

%% Conservative
for i=1:N_sim
    L1_cons(i) = sum(abs(N_cons_all(i,:)- N_ANA_all(i,:)))/sum(abs(N_ANA_all(i,:)));
    L2_cons(i) = sum((N_cons_all(i,:)- N_ANA_all(i,:)).^2)/sum((N_ANA_all(i,:)).^2);
    if i~=N_sim
        EOC_cons(i) = log(sum(abs(N_cons_all(i,:)- N_ANA_all(i,:)))/sum(abs(N_cons_all(i+1,:)- N_ANA_all(i+1,:))))/log(2);
    end
end

%% NON LINEAR BREAKAGE CODE (2-6-21) (ASHOK DAS) - 2D
 clear all
 close all

example = 4;
%% Inputs
x1_min = 0; x1_max = 1; I1 = 20; % 1st property limits
x2_min = 0; x2_max = 1; I2 = 20; % 2nd property limits
x3_min = 0; x3_max = 1; I3 = 20; % 3rd property limits

T =10;  % [sec] Process time
len_T = 6;
time  = linspace(0,T,len_T); % Time discretization

K_index = 2; % 1-> K=1; 2-> K=x1*x2*y1*y2

%% Initial PSD
N_ini    = zeros(I1*I2*I3,1); % Initialization
N_ini(I1*I2*I3) = 1;

tic
%% Discretizations and other functions
[x1,R1,del_x1] = Grids2(x1_min, x1_max, I1);
[x2,R2,del_x2] = Grids2(x2_min, x2_max, I2);
[x3,R3,del_x3] = Grids2(x3_min, x3_max, I3);% x-> pivot pts; R-> boundary pts; del_x-> grid length
% [x1,R1,del_x1] = Lin_Grids(x1_min, x1_max, I1);
% [x2,R2,del_x2] = Lin_Grids(x2_min, x2_max, I2); % x-> pivot pts; R-> boundary pts; del_x-> grid length
fprintf('Discretization done...')
%% K function
K = zeros(I1,I2,I3,I1,I2,I3); % initialization
i=1:I1; j=1:I2;

for c1=1:I3
    for c2=1:I1
        for c3=1:I2
            for c4=1:I3
                K(:,:,c1,c2,c3,c4) = (x1(i)'*x2(j))*(x3(c1)*x1(c2)*x2(c3)*x3(c4));
            end
        end
    end
end
fprintf('K-funtion done...')

%% Beta function
p1 = p_Fun_mat(x1,R1,I1); % p(i,m) matrix form (1st variable)
p2 = p_Fun_mat(x2,R2,I2); % p(i,m) matrix form (2nd variable)
p3 = p_Fun_mat(x3,R3,I3); % p(i,m) matrix form (3rd variable)


B = B_Fun(p1,p2,p3,x1,x2,x3,R1,R2,R3); % matrix version of \int_{R1(i1)}^{p1(i1,m1)} \int_{R2(i2)}^{p2(i2,m2)} b(x1,x2;x1_m1,x2_m2) dx1 dx2; where b=4/x1(m1)*x2(m2)
fprintf('Beta function done...')

%% breakage distribution function  (FOR CONSERVATIVE APPROACH-1: MID-PT. RULE)
% beta = zeros(I1,I2,I3,I1,I2,I3); % Initialization
% for i=1:I1
%     for j=1:I2
%         for k=1:I3
%             beta(:,:,:,i,j,k) = 8/(x1(i)*x2(j)*x3(k));
%         end
%     end
% end

% Int. of b(x,y;z) = 4/y1*y2   (FOR CONSERVATIVE APPROACH-2: ANALYTIC INT.)
% beta_cons = zeros(I1,I2,I1,I2); % Initialization
% for i=1:I1
%     for j=1:I2
%         beta_cons(:,:,i,j) = 4*log(R1(i+1)/R1(i))*log(R2(j+1)/R2(j)) /(del_x1(i)*del_x2(j));
%     end
% end
fprintf('br. dist function done...')

%% Functions related to schemes
[w1,w2_b,w2_d] = weights(x1,x2,x3,B); % Weight functions for MC and NPMC
fprintf('Weight creation done...')

t_sim_0 = toc;
%% Solution part
options = odeset('RelTol',1e-6, 'AbsTol',1e-6);
%options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'MaxSTep',1000,'Stat','on');
tic
[T1,N1] = ode45(@discrete_MC, time, N_ini, options, K,B,w1,x1,x2,x3); % Mass conserving technique
t_sim_1 = toc;

tic
[T2,N2] = ode45(@discrete_NPMC, time, N_ini, options, K,B,w2_b,w2_d,x1,x2,x3); % Number conserve + Mass conserving technique
t_sim_2 = toc;

% tic
% [T3,N3] = ode45(@discrete_conserve, time, N_ini, options, x1,x2,del_x1,del_x2,K,beta); % Conservative approach-1
% t_sim_3 = toc;
% 
% tic
% [T4,N4] = ode45(@discrete_conserve, time, N_ini, options, x1,x2,del_x1,del_x2,K,beta_cons); % Conservative approach-2
% t_sim_4 = toc;


%% Final PSD
% PSD_MC = vec2mat(N1,I1,I2);  PSD_NPMC = vec2mat(N2,I1,I2); PSD_cons = vec2mat(N4,I1,I2);
% [x,y]=meshgrid(x1,x2);
%
% figure
% surf(log(x),log(y),PSD_MC)
% figure
% surf(x,y,PSD_NPMC)
% figure
% surf(x,y,PSD_cons)

%% Figure plot
%% Total particle - M_00
N_tot1 = sum(N1,2); N_tot2 = sum(N2,2); %N_tot3 = sum(N3,2); N_tot4 = sum(N4,2);
figure
plot(T1,N_tot1,'bo--','linewidth',1.5,'markersize',11)
hold on
plot(T2,N_tot2,'rs--','linewidth',1.5,'markersize',11)
% plot(T3,N_tot3,'+-','markersize',12)
%plot(T4,N_tot4,'m^--','linewidth',1.5,'markersize',11)
% Analytical result
plot(time,(7*time +1),'k-','markersize',16,'linewidth',2.5)

legend({'HC','HCNP','Exact'},'fontsize',18,'Location','northwest') %'CFV',
xlabel('Time (t)','fontsize',25);
ylabel('M_{0,0,0}(t)','fontsize',25);
% savePDF(['Ex_',num2str(example),'_M000'])

%% Total mass - M_{1,1,1}
x1x2 = x1'*x2; Hyp_vol_ind_mat = zeros(I1,I2,I3);
for k=1:I3
    Hyp_vol_ind_mat(:,:,k) = x1x2 * x3(k);
end

Hyp_vol_ind = reshape(Hyp_vol_ind_mat,[],1);

M_tot_1 = N1*Hyp_vol_ind; M_tot_2 = N2*Hyp_vol_ind;
%M_tot_3 = N3*area_mat_vec; M_tot_4 = N4*area_mat_vec;

figure
plot(T1,M_tot_1,'bo','linewidth',1.5,'markersize',18)
hold on
plot(T2,M_tot_2,'rs','linewidth',1.5,'markersize',14)
% plot(T3,M_tot_3,'+-','markersize',12)
%plot(T4,M_tot_4,'m^','linewidth',1.5,'markersize',23)

plot(time,ones(len_T,1),'k-','markersize',16,'linewidth',2.5)
% title('Total mass')
legend({'HC','HCNP','Exact'},'fontsize',18,'Location','best') %'CFV',
ylim([0.998 1.003])
xlabel('Time (t)','fontsize',25);
ylabel('M_{1,1,1}(t)','fontsize',25);
% savePDF(['Ex_',num2str(example),'_M111'])

%% M_{1,0,0} + M_{0,1,0}+ M_{0,0,1}
Mix_1 = zeros(1,len_T); Mix_2 = zeros(1,len_T);
for l=1:len_T
    N_MC = reshape(N1(l,:),I1,I2,I3); N_NPMC = reshape(N2(l,:),I1,I2,I3);
    for i1=1:I1
        for i2=1:I2
            for i3=1:I3
                Mix_1(l) = Mix_1(l)+ (x1(i1)+x2(i2)+x3(i3))*N_MC(i1,i2,i3);
                Mix_2(l) = Mix_2(l)+ (x1(i1)+x2(i2)+x3(i3))*N_NPMC(i1,i2,i3);
            end
        end
    end
end


% for i=1:len_T
%     N_MC = vec2mat(N1(i,:),I1,I2); N_NPMC = vec2mat(N2(i,:),I1,I2);
%     N_cons1 = vec2mat(N3(i,:),I1,I2); N_cons2 = vec2mat(N4(i,:),I1,I2); 
%     
%     Mix_1(i) = x1*sum(N_MC,2) + sum(N_MC)*x2'; Mix_2(i) = x1*sum(N_NPMC,2) + sum(N_NPMC)*x2';
%     Mix_3(i) = x1*sum(N_cons1,2) + sum(N_cons1)*x2'; Mix_4(i) = x1*sum(N_cons2,2) + sum(N_cons2)*x2';
% end
figure
plot(T1,Mix_1,'bo--','linewidth',1.5,'markersize',11)
hold on
plot(T2,Mix_2,'rs--','linewidth',1.5,'markersize',11)
% plot(T3,Mix_3,'+-','markersize',12)
%plot(T4,Mix_4,'m^--','linewidth',1.5,'markersize',11)

legend({'HC','HCNP','CFV'},'fontsize',18,'Location','northwest')
%ylim([0.998 1.003])
xlabel('Time (t)','fontsize',25);
ylabel('M_{1,0,0}(t) + M_{0,1,0}(t) +  M_{0,1,0}(t)','fontsize',22);
% savePDF(['Ex_',num2str(example),'_M100_M010_M001'])

%% Average hypervolume
figure
plot(T1,M_tot_1./N_tot1,'bo--','linewidth',1.5,'markersize',11)
hold on
plot(T2,M_tot_2./N_tot2,'rs--','linewidth',1.5,'markersize',11)
% plot(T3,M_tot_3./N_tot3,'+-','markersize',12)
%plot(T4,M_tot_4./N_tot4,'m^--','linewidth',1.5,'markersize',11)

plot(time,1./(7*time+1),'k-','markersize',16,'linewidth',2.5)
% title('Total mass')
legend({'HC','HCNP','Exact'},'fontsize',18,'Location','best') % ,'CFV'
xlabel('Time (t)','fontsize',25);
ylabel('Average hypervolume','fontsize',25);
% savePDF(['Ex_',num2str(example),'_avg_hypervol'])


%% Error table

M000_err_1 = ((7*time+1) - N_tot1')./ (7*time+1); M000_err_2 = ((7*time+1) - N_tot2')./ (7*time+1);
%M00_err_3 = ((3*time+1) - N_tot3')./ (3*time+1); M00_err_4 = ((3*time+1) - N_tot4')./ (3*time+1);

M111_err_1 = (1 - M_tot_1'); M111_err_2 = (1 - M_tot_2');
%M111_err_3 = (1 - M_tot_3'); M11_err_4 = (1 - M_tot_4');

% Mix_err_1 = (2 - Mix_1)/2; Mix_err_2 = (2 - Mix_2)/2;
% Mix_err_3 = (2 - Mix_3)/2; Mix_err_4 = (2 - Mix_4)/2;
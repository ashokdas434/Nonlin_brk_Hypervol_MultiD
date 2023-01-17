%% NON LINEAR BREAKAGE CODE (2-6-21) (ASHOK DAS) - 1D
%clear all
%close all
global nu
nu = 1;
example = 1; % 1 or 2   % example number

%% ***************************** Inputs *************************
for z= 0%:5
    x_min = 0; x_max = 1;  % property limits
    I     = 30; % Number of intervals in base case
    grid_mul = 2^z; % Power of 2

    if example == 1
        T = 100;  % [sec] Process time
    else
        T = 1.75;  %0.4;%   [sec] Process time
    end
    len_T = 101;
    time  = linspace(0,T,len_T); % Time discretization

    %% Discretizations and other functions
    [x,R,del_x,L] = Grids2(x_min, x_max, I, grid_mul); % x-> pivot pts; R-> boundary pts; del_x-> grid length
    %[x,R,del_x] = Lin_Grids(x_min, x_max, I); % x-> pivot pts; R-> boundary pts; del_x-> grid length

    p = p_Fun_mat(x,R,L); % p(i,m) matrix form

    %% Initial PSD
    N_ini    = zeros(L,1);  N_ini(L) = 1;  % Initialization

    %% Different functions for each example
    K = K_Fun(example,x,L); % K-function (collision freq function)
    [b_Fun, b_mat, beta, beta_cons, frag] = b_Function(example,x,p,R); % breakage dist fun related matrix and fun_handles

    %% Functions related to schemes
    [w1,w2_b,w2_d] = weights(x,beta, frag); % Weight functions for MC and NPMC

    %% **************** Numerical Solution part *****************************
    options = odeset('RelTol',1e-6, 'AbsTol',1e-6);
    %options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'MaxSTep',1000,'Stat','on');
    tic
    [T1,N1] = ode45(@discrete_MC, time, N_ini, options, K,beta,w1,x); % Mass conserving technique
    t_sim_1 = toc;

    tic
    [T2,N2] = ode45(@discrete_NPMC, time, N_ini, options, K,beta,w2_b,w2_d,x); % Number conserve + Mass conserving technique
    t_sim_2 = toc;

    tic
    [T3,N3] = ode45(@discrete_conserve, time, N_ini, options, x,del_x,K,b_mat); % Conservative approach-1
    t_sim_3 = toc;

    tic
    [T4,N4] = ode45(@discrete_conserve, time, N_ini, options, x,del_x,K,beta_cons); % Conservative approach-2
    t_sim_4 = toc;

    %% ********************** Analytical solution *******************************
    [N_ana, N_tot_ana, M_dist_ana, M_tot_ana, M2_tot_ana] = Analytical_sol(example,x,R,time);

    %% ************************ data saving *********************************
    N_MC = N1(end,:); N_NPMC = N2(end,:); N_cons = N4(end,:); N_ANA = N_ana(end,:);
   % save(['sim2_Ex-',num2str(example),'-Grid_mul-',num2str(grid_mul),'.mat'],'x','del_x','N_MC','N_NPMC','N_cons','N_ANA')



    %% *************************** Figure plot ***********************************
    %% PSD
    %for i = len_T-20:4:len_T
    i=len_T;
    figure
    semilogx(x,N1(i,:),'bo','linewidth',2.5,'markersize',11)
    hold on
    semilogx(x,N2(i,:),'rs','linewidth',2.5,'markersize',11)
%    semilogx(x,N3(i,:),'^-','markersize',12)
    semilogx(x,N4(i,:),'m^','linewidth',2.5,'markersize',11)

    semilogx(x,N_ana(i,:),'k-','markersize',16,'linewidth',2.5)
  %  title('PSD')
    legend({'HC','HCNP','CFV','Exact'},'fontsize',18,'Location','best')
    xlim([0 1])
    xlabel('Particle size','fontsize',25);
    ylabel('Number of fragments','fontsize',25);
    savePDF(['Ex_',num2str(example),'_PSD'])
    % end

    %% PSD (number density function)
%     %for i = len_T-20:4:len_T
%     i=len_T;
%     figure
%     semilogx(x,N1(i,:)./del_x,'o-')
%     hold on
%     semilogx(x,N2(i,:)./del_x,'*-')
%     semilogx(x,N3(i,:)./del_x,'^-','markersize',12)
%     semilogx(x,N4(i,:)./del_x,'p-','markersize',12)
% 
%     semilogx(x,N_ana(i,:)./del_x,'-','markersize',16,'linewidth',2)
%     title('PSD')
%     legend('MC','NPMC','Conserve','Conserve 2','Analytical')
%     % end

    %% Total particle
    N_tot1 = sum(N1,2); N_tot2 = sum(N2,2); N_tot3 = sum(N3,2); N_tot4 = sum(N4,2);
    figure
    plot(T1(1:10:end),N_tot1(1:10:end),'bo','linewidth',2.5,'markersize',13)
    hold on
    plot(T2(1:10:end),N_tot2(1:10:end),'rs','linewidth',2.5,'markersize',13)
 %   plot(T3,sum(N3,2),'+-','markersize',12)
    plot(T4(1:10:end),N_tot4(1:10:end),'m^','linewidth',2.5,'markersize',13)

    plot(time(1:10:end),N_tot_ana(1:10:end),'k-','markersize',16,'linewidth',2.5)
    %title('No of partile')
    legend({'HC','HCNP','CFV','Exact'},'fontsize',18,'Location','northwest')
%     xlim([0 1])
    xlabel('Time (t)','fontsize',25);
    ylabel('M_0(t)','fontsize',25);
    savePDF(['Ex_',num2str(example),'_M0'])

    %% Total Mass
    M_tot1 = N1*x'; M_tot2 = N2*x'; M_tot3 = N3*x'; M_tot4 = N4*x';
    figure
    plot(T1(1:10:end),M_tot1(1:10:end),'bo','linewidth',1.5,'markersize',18)
    hold on
    plot(T2(1:10:end),M_tot2(1:10:end),'rs','linewidth',1.5,'markersize',14)
%    plot(T3,N3*x','+-','markersize',12)
    plot(T4(1:10:end),M_tot4(1:10:end),'m^','linewidth',1.5,'markersize',23)

    plot(time,M_tot_ana,'k-','markersize',16,'linewidth',2.5)
   % title('Total mass')
    legend({'HC','HCNP','CFV','Exact'},'fontsize',18)
     ylim([0.998 1.003])
    xlabel('Time (t)','fontsize',25);
    ylabel('M_1(t)','fontsize',25);
    savePDF(['Ex_',num2str(example),'_M1'])

      %% M_2
      xx = x.*x;
      M2_tot1 = N1*xx'; M2_tot2 = N2*xx'; M2_tot3 = N3*xx'; M2_tot4 = N4*xx';
      l = 0;
    figure
    plot(T1([1:4 6:5:41]),M2_tot1([1:4 6:5:41]),'bo','linewidth',2.5,'markersize',11)
    hold on
    plot(T2([1:4 6:5:41]),M2_tot2([1:4 6:5:41]),'rs','linewidth',2.5,'markersize',11)
%    plot(T3,N3*xx','+-','markersize',12)
    plot(T4([1:4 6:5:41]),M2_tot4([1:4 6:5:41]),'m^','linewidth',2.5,'markersize',11)

    plot(time([1:4 6:5:41]),M2_tot_ana([1:4 6:5:41]),'k-','markersize',16,'linewidth',2.5)

    legend({'HC','HCNP','CFV','Exact'},'fontsize',18)
%     xlim([0 1])
    xlabel('Time (t)','fontsize',25);
    ylabel('M_2(t)','fontsize',25);
    savePDF(['Ex_',num2str(example),'_M2'])

end
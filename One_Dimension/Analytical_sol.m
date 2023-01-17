%% Function to calculate analytical solutions
function [N_ana, N_tot_ana, M_dist_ana, M_tot_ana, M2_tot_ana] = Analytical_sol(example,x,R,time)

global nu

I = length(x); len_T = length(time);
N_ana = zeros(len_T,I); M_dist_ana = zeros(len_T,I); % Initialize

switch example
    case 1
        f = @(x,t) exp(-x*t).*(2.*t+ t.^2 .*(1-x));
        xf = @(x,t) exp(-x*t).*(2.*t.*x+ t.^2 .*(x-x.^2));
        xxf = @(x,t) exp(-x*t).*(2.*t.*x.^2+ t.^2 .*(x.^2-x.^3));

        for j=1:len_T
            for i=1:I-1
                N_ana(j,i) = integral(@(x) f(x,time(j)),R(i),R(i+1));
                M_dist_ana(j,i) = integral(@(x) xf(x,time(j)),R(i),R(i+1));
                M2_dist_ana(j,i) = integral(@(x) xxf(x,time(j)),R(i),R(i+1));

            end
            N_ana(j,I) = integral(@(x) f(x,time(j)),R(I),x(I)) + exp(-time(j));
            M_dist_ana(j,I) = integral(@(x) xf(x,time(j)),R(I),x(I)) + exp(-time(j));
            M2_dist_ana(j,I) = integral(@(x) xxf(x,time(j)),R(I),x(I)) + exp(-time(j));
        end
        
        
        
    case 2
        b0 = (nu +2)/(nu+1); 
        tau = @(t) log(1-(b0-1)*t)/(1-b0);
        f = @(x,t) exp(-tau(t)).* x.^nu .* (tau(t).*(nu+2)./log(1./x)).^(1/2) .*besseli(1,2*(tau(t).*(nu+2).*log(1./x)).^(1/2));
        xf = @(x,t) exp(-tau(t)).* x.^(nu+1) .* (tau(t).*(nu+2)./log(1./x)).^(1/2) .*besseli(1,2*(tau(t).*(nu+2).*log(1./x)).^(1/2));
        xxf = @(x,t) exp(-tau(t)).* x.^(nu+2) .* (tau(t).*(nu+2)./log(1./x)).^(1/2) .*besseli(1,2*(tau(t).*(nu+2).*log(1./x)).^(1/2));
        
        for j=1:len_T
            for i=1:I-1
                N_ana(j,i) = integral(@(x) f(x,time(j)),R(i),R(i+1));
                M_dist_ana(j,i) = integral(@(x) xf(x,time(j)),R(i),R(i+1));
                M2_dist_ana(j,i) = integral(@(x) xxf(x,time(j)),R(i),R(i+1));
            end
            N_ana(j,I) = integral(@(x) f(x,time(j)),R(I),x(I)) + exp(-tau(time(j)));
            M_dist_ana(j,I) = integral(@(x) xf(x,time(j)),R(I),x(I)) + exp(-tau(time(j)));
            M2_dist_ana(j,I) = integral(@(x) xxf(x,time(j)),R(I),x(I)) + exp(-tau(time(j)));
        end
end


N_tot_ana = sum(N_ana,2); % Total particles
M_tot_ana = sum(M_dist_ana,2); % Total mass
M2_tot_ana = sum(M2_dist_ana,2); % Total mass
return
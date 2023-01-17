%% breakage dist fun in function handle and matrix form
function [b_Fun, b_mat, beta, beta_cons, frag] = b_Function(example,x,p,R)

global nu
I = length(x);
b_mat = zeros(I); beta = zeros(I); beta_cons = zeros(I); % initialization

switch example
    case 1 % b(x,y) = 2/y
        b_Fun = @(a,b) ones(1,length(a)).*2./b; % for 1D int
        b_Fun_cons = @(a,b) 2./b; % for conserve form (2d int)
        
        for i=1:I
            for j=1:I
                b_mat(i,j) = b_Fun(x(i),x(j));
                beta(i,j) = integral(@(a) b_Fun(a,x(j)),R(i),p(i,j));
                
                del = (R(i+1)-R(i))*(R(j+1)-R(j));
                beta_cons(i,j) = integral2(b_Fun_cons,R(i),R(i+1),R(j),R(j+1)) /del;
            end
        end
        
    case 2 % K_st(i,j) = i*j
        b_Fun = @(a,b) (nu+2).*a.^(nu)./b.^(1+nu);
        
        for i=1:I
            for j=1:I
                b_mat(i,j) = b_Fun(x(i),x(j));
                beta(i,j) = integral(@(a) b_Fun(a,x(j)),R(i),p(i,j));
                
                del = (R(i+1)-R(i))*(R(j+1)-R(j));
                beta_cons(i,j) = integral2(b_Fun,R(i),R(i+1),R(j),R(j+1)) /del;
            end
        end
        

end

%% No. of fragments per breakage
Y=10; % random value
frag = integral(@(x) b_Fun(x,Y),0,Y);

return
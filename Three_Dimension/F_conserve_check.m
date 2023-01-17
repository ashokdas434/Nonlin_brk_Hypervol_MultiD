%% Function to calculate F for consevative approach
function [F1,F2,F3] = F_conserve_check(x1,x2,del_x1,del_x2,N_mat,K,beta)

I1 = length(x1); I2 = length(x2);
F1 = zeros(I1+1,I2+1); F2 = zeros(I1+1,I2+1); F3 = zeros(I1+1,I2+1);% Initialization

x1_del_x1 = x1.*del_x1; x2_del_x2 = x2.*del_x2;
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
        s1=0;   %% F1
        for m2=i+1:I1
            for n2=j:I2
                s1 = s1+ NKN(m2,n2)*(x1_del_x1(1:i)*beta(1:i,j,m2,n2));
            end
        end
        F1(i+1,j)= x2(j)*s1;
        
        s2=0;  %%F2
        for m1=i:I1
            for m2=(j+1):I2
                s2 = s2+ NKN(m1,m2)*(x2_del_x2(1:j)*beta(i,1:j,m1,m2)');
            end
        end
        F2(i,j+1)= x1(i)*s2;
        
        s3=0;  %%F3
        for m1=(i+1):I1
            for m2=(j+1):I2
                temp = beta(1:i,1:j,m1,m2)*x2_del_x2(1:j)';                
                s3 = s3+ NKN(m1,m2)*(x1_del_x1(1:i)*temp);
            end
        end
        F3(i+1,j+1)= s3;
    end
end
        

return

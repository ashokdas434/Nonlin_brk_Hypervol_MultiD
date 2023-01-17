%% Function to formulate the weights
function [w1,w2_b,w2_d] = weights(x1,x2,x3,B)

I1 = length(x1); I2 = length(x2); I3 = length(x3);
w1 = zeros(I1,I2,I3); w2_b = zeros(I1,I2,I3); w2_d = zeros(I1,I2,I3); % Initialization
nu = 8; % number of fragments per breakage

x1x2 = x1'*x2;  x1x2x3= zeros(I1,I2,I3);
for c2=1:I3
    x1x2x3(:,:,c2) = x3(c2) * x1x2;
end

for i1=1:I1
    for i2=1:I2
        for i3=1:I3
            j1=1:i1;  j2=1:i2;  j3=1:i3;
            S=0;
            for c1=1:i3
                S = S + x3(c1) * (x1(j1)* (B(j1,j2,c1,i1,i2,i3)*x2(j2)') );
            end
            w1(i1,i2,i3) = S/(x1(i1)*x2(i2)*x3(i3));

            temp1 = x1(i1)*x2(i2)*x3(i3)- x1x2x3(j1,j2,j3);
            temp2 = temp1.* B(j1,j2,j3,i1,i2,i3);
            S2 = sum(sum(sum(temp2)));
            w2_b(i1,i2,i3) = x1(i1)*x2(i2)*x3(i3)* (nu -1)/S2;
            w2_d(i1,i2,i3) = w2_b(i1,i2,i3) * w1(i1,i2,i3);
        end
    end
end

w2_b(1,1,1) = 0; w2_d(1,1,1) = 0;

return


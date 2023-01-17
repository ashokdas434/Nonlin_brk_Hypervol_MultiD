%% This function creates the matrix version of
% \int_{R1(i1)}^{p1(i1,m1)} \int_{R2(i2)}^{p2(i2,m2)} b(x1,x2;x1_m1,x2_m2) dx1 dx2; where b=4/x1(m1)*x2(m2)
function B = B_Fun(p1,p2,p3,x1,x2,x3,R1,R2,R3)
I1 = length(x1); I2 = length(x2); I3 = length(x3);
B = zeros(I1,I2,I3,I1,I2,I3); % initialization

for i1=1:I1
    for i2=1:I2
        for i3=1:I3
            for m1=1:I1
                for m2=1:I2
                    for m3=1:I3
                        B(i1,i2,i3,m1,m2,m3) = 8*(p1(i1,m1)-R1(i1))*(p2(i2,m2)-R2(i2))*(p3(i3,m3)-R3(i3))/(x1(m1)*x2(m2)*x3(m3));
                    end
                end
            end
        end
    end
end

return
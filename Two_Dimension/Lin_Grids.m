function [x,R,del_x] = Lin_Grids(x_min, x_max, I)

del_x = (x_max - x_min)/I * ones(1,I);

R = linspace(x_min,x_max,I+1);
R(I+1) = R(I+1) + (R(I+1) -R(I));

x = zeros(1,I);
for i=1:I
    x(i) = 0.5*(R(i)+R(i+1));
end

end
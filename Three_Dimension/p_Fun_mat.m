%% This function determines the p-value (matrix form)
function p = p_Fun_mat(x,R,I)
p = zeros(I); % initialization

for i=1:I
    p(i,:) = R(i+1);
end

for i=1:I
    p(i,i) = x(i);
end
  
return
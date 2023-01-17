%% This function determines the p-value
function res = p_Fun(i,m,x,R)

if i==m
    res = x(i);
else
    res = R(i+1);
end

return
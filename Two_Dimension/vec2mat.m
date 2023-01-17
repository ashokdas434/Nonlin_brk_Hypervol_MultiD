% Function: Vector -> Matrix
function mat = vec2mat(vec,I1,I2)

mat = zeros(I1,I2);

for i=1:I1
    mat(i,:) = vec((i-1)*I2+1:i*I2);
end

return
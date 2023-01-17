% Function: Matrix -> Vector
function vec = mat2vec(mat)

[I1,I2] = size(mat);
vec = zeros(I1*I2,1);

for i=1:I1
    vec((i-1)*I2+1:i*I2) = mat(i,:);
end

return
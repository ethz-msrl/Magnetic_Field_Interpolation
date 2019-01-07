function [ sum ] = tricubic( A, x, y, z )
%TRICUBIC Summary of this function goes here
%   Detailed explanation goes here

sum = 0;
for i = 0:3
    for j=0:3
        for k=0:3
            sum = sum + A(i+1,j+1,k+1) * x^i * y^j * z^k;
        end
    end
end

end
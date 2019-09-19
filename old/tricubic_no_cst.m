function [ sum ] = tricubic_no_cst( av, x, y, z )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
sum = 0;
ind = 1;
for i = 0:3
    for j=0:3
        for k=0:3
            if (i == 0) && (j == 0) && (k==0)
                continue
            end

            sum = sum + av(ind) * x^i * y^j * z^k;
            
            ind = ind + 1;
        end
    end
end

end


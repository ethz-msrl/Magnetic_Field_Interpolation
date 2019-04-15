function [ y ] = tricubic_grad( av, x, y, z )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
sumx = 0;
sumy = 0;
sumz = 0;

ind = 1;
for i = 0:3
    for j=0:3
        for k=0:3
            if (i == 0) && (j == 0) && (k==0)
                continue
            end

            if (i > 0)
                sumx = sumx + i * av(ind) * x^(i-1) * y^j * z^k;
            end
            
            if (j > 0)
                sumy = sumy + j * av(ind) * x^i * y^(j-1) * z^k;
            end
            if (k > 0)
                sumz = sumz + k * av(ind) * x^i * y^j * z^(k-1);
            end
            
            ind = ind + 1;
        end
    end
end

y = [sumx; sumy; sumz];

end


function [b,a] = cantor_pairing_1to2(z)
w = floor((sqrt(8.*z+1)-1)./2);
t = (w.^2 + w)./2;
a = z-t;
b = w-a;
end
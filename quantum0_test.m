clear all;
p = [0;0;0.5;0;0;0;-0.5i;0;-0.5;0;0;0;0;0;0.5i;0]; 

f = ones(16,16);
w = exp((i*pi/8));

for i = 1:16
    for j = 1:16
        f(i,j) = w^((i-1)*(j-1));
    end
end

a = (f*p)/4;
x = exp(i*pi/8);
t=[];
for k = 1:100
    x1 = x^k;
    t=[t;x1];
end
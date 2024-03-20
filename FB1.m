function [a3] = FB1(f,alfa)
N=length(f);
nb=(1:N);
for m1=1:N
a3(m1)=(2/(N^2*(besselj(1,alfa(m1))).^2))*sum(nb.*f.*besselj(0,alfa(m1)/N*nb));
end
function [alfa] = root1(f)
N=length(f);
%computation of roots of bessel function Jo(x)
if exist('alfa') == 0
    x=2;
    alfa=zeros(1,N);
    for i=1:N
        ex=1;
        while abs(ex)>.00001
            ex=-besselj(0,x)/besselj(1,x);
            x=x-ex;
        end
        alfa(i)=x;
        x=x+pi;
    end
end
end
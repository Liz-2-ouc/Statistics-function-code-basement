function [s, T, s_alf ] = power_discrt(x,alpha)
%输入参数x为一维时间序列， 
%输出参数s, T 和s_alf为长度相同的一维数组，分别代表不同频率的功率谱值、周期、和临界谱值
n = length(x);
a=zeros(n);b=zeros(n);
for k = 1:1:(n-1)/2
    a1=0;b1=0;
    for t = 1:1:n
        a1 = a1+x(t)*cos(2*pi*k*t/n);
        b1 = b1+x(t)*sin(2*pi*k*t/n);
        a2 = x(t)*cos(pi*t);
    end
    a(k) = a1*2/n;b(k) = b1*2/n;
    s(k) = (a(k)^2+b(k)^2)/2;
    T(k) = n/k;
end
    if mod(n,2)==0
        a(n/2) = a2/n;
        b(n/2) = 0;
    end
%%显著性检验
Fc=finv(1-alpha/2,2,n-2-1);
sc=2*Fc*var(x)/(2*Fc+(n-2-1));
s_alf(1:length(T))=sc;
end
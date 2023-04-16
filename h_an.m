function [a_0,a_k,b_k,phi,A_k,S_k,vf]=h_an(x,type)
%本函数用于傅氏系数的求解计算
%A_k为各谐波振幅；a_0为第零傅氏系数;a_k为第k级余弦傅氏系数;b_k为第k级正弦傅氏系数;
%phi为各谐波初相位;vf为各方差贡献率;S_k为各谐波功率 

ss=std(x,1).^2;% ss为时间序列方差
%计算谐波数P
n=length(x);
if mod(n,2)==1
    p=(n-1)./2;
else
    p=n./2;
end
%初始化各系数
a_0=mean(x);
a_k=zeros(p,1);
b_k=zeros(p,1);
phi=zeros(p,1);
A_k=zeros(p,1);
S_k=zeros(p,1);
for i=1:p
    if type==3
        %使用第三种形式时：
        for j=1:n
            a_k(i)=a_k(i)+2./n.*(x(j).*cos(i.*2.*pi./n.*j));
            b_k(i)=b_k(i)+2./n.*(x(j).*sin(i.*2.*pi./n.*j));
        end
        
        if mod(n,2)==0
            %p为偶数时不具有双边性
            a_k(p)=0;%初始化
            for j=1:n
                a_k(p)=a_k(p)+1./n.*(x(j).*cos(pi.*j));
                b_k(p)=0;
            end
        end
    elseif type==2
        %使用第二种形式时：
        for j=1:n
            a_k(i)=a_k(i)+2./n.*(x(j).*cos(i.*2.*pi./n.*(j-1)));
            b_k(i)=b_k(i)+2./n.*(x(j).*sin(i.*2.*pi./n.*(j-1)));
        end
        
        if mod(n,2)==0
            %p为偶数时不具有双边性
            a_k(p)=0;%初始化
            for j=1:n
                a_k(p)=a_k(p)+1./n.*(x(j).*cos(pi.*(j-1)));
                b_k(p)=0;
            end
        end
    else
        %使用第一种形式时：
        for j=1:n
            a_k(i)=a_k(i)+2./n.*(x(j).*cos(i.*2.*pi./n.*(j-1)));
            b_k(i)=b_k(i)+2./n.*(x(j).*sin(i.*2.*pi./n.*(j-1)));
        end
        
        if mod(n,2)==0
            %p为偶数时不具有双边性
            a_k(p)=0;%初始化
            for j=1:n
                a_k(p)=a_k(p)+1./n.*(x(j).*cos(pi.*(j-1)));
                b_k(p)=0;
            end
        end
    end
end

for i=1:p
    phi(i) = angle(complex(b_k(i),a_k(i)));
    A_k(i)=sqrt(a_k(i).^2+b_k(i).^2);
    S_k(i)=1./2.*A_k(i).^2;
    vf(i)=1./2.*A_k(i).^2./ss;
end

if mod(n,2)==0
    S_k(p)=A_k(p).^2;
    vf(i)=A_k(p).^2./ss;
end
   
end
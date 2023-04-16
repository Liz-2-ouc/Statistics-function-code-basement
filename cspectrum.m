function [T,S_l,strw,strw_a]=cspectrum(x,m,a,sort);
%m代表最大滞后相关长度，x代表原始数据，a代表检验所需要的置信度
%sort 代表排序方式(均为从大到小)，1：波数；2：频率；3：周期
%S_l代表连续功率谱，strw代表红白噪声功率谱，strw_a代表a置信度的红白噪声上限

n=length(x);%时间序列长度

%% ------------对时间序列进行标准化处理-------------
x_ave=mean(x);
x_s=std(x,1);
x_z=(x-x_ave)./x_s;

%% ------------计算自相关函数-------------
r=zeros(1,m);
% 其中t代表时滞
for t=1:m
    r(t)=dot(x_z(1:end-t),x_z(1+t:end))./(length(x_z(1:end-t)));
end
r0=x_z*x_z'./length(x_z);

%% ------------功率谱计算-------------
S_l=zeros(m+1,1);
%l代表波数
for l=0:m
    ss=0;
    % 其中t代表时滞
    for t=1:(m-1)
        ss=ss+r(t).*(1+cos(pi.*t./m)).*cos(l.*pi.*t./m); %这一步已经进行了平滑处理
    end
    
    if l==0 || l==m
        B=0.5;
    else 
        B=1;
    end
    S_l(l+1,1)=(r0+ss)*B./m;
end

%% ------------功率谱计算-------------
%根据t分布的上分位点表得到判据Ra
Ra=(-1+tinv(1-a,n-2).*sqrt(n-2))./(n-1);
if Ra(1) <= r(1)
    disp('应用红噪声谱检验') 
    w=zeros(m+1,1);%不同周期对应的波长频率 
    strw=zeros(m+1,1);strw_a=zeros(m+1,1);
    for l=0:m
        j=l+1;
        dt=1;%设置时间间隔为1
        w(j)=pi*l*dt/m;     
        strw(j)=mean(S_l).*(1-r(1)^2)/(1-2.*r(1).*cos(w(j))+r(1)^2);  
        z=(2*n-m/2)/m;
        X_2=chi2inv(1-a,z);  
        strw_a(j)=X_2.*strw(j)./z;  
    end
else
    disp('应用白噪声谱检验')   
    w=zeros(m+1,1);%不同周期对应的波长频率 
    strw=zeros(m+1,1);strw_a=zeros(m+1,1);
    for l=0:m
        j=l+1;
        dt=1;%设置时间间隔为1
        w(j)=pi*l*dt/m;     
        strw(j)=mean(S_l); 
        z=(2*n-m/2)/m;
        X_2=chi2inv(1-a,z);  
        strw_a(j)=X_2.*strw(j)./z;    
    end
end

l=[0:m];
T=2*m./l;

if sort==3
     T=flip(T);
     S_l=flip(S_l);
     strw=flip(strw);
     strw_a=flip(strw_a);
end   

end
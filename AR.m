function phi=AR(x,p)
    %标准化
    x_ave=mean(x);
    x_s=std(x,1);
    x_z=(x-x_ave)./x_s;
    %得到自相关系数函数，其中滞后系数从1到p
    %先进行初始化
    r=zeros(1,p);
    for t=1:length(x)-1
        r(t)=dot(x_z(1:end-t),x_z(1+t:end))./(length(x_z(1:end-t)));
    end
    %建立Yule-Walker方程并求解自回归系数
    S=ones(p);
    for i=1:p-1
        for j=1:p-1
            S(i+j,i)=r(j);
            S(i,i+j)=r(j);
        end
    end
    S=S(1:p,1:p);
    b=r(1:p);
    phi=S\b';
end
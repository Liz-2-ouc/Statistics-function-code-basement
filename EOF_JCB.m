function [L,Y,lambda,vf,cvf]= EOF_JCB(X,ks,kks)
%输入变量说明
%X为原资料矩阵[m n]，m为变量数，n为样本长度
%ks为资料变量设置，ks = -1采用原始阵；ks = 0中心化阵；ks = 1采用标准化矩阵
%kks为是否赋物理意义选项，kks = 1表示赋物理意义；kks = 0表示不赋物理意义
%
%输出变量说明
%L为系数矩阵[m mnl]，Y为主成分矩阵[mnl n]，lambda为特征向量
%vf为方差贡献率，cvf为累积方差贡献率

%% 判断资料变量设置
if ks == 0% 距平阵
    X = X - mean(X,2);
elseif ks == 1% 标准化矩阵
    X =  (X - mean(X,2))./std(X,0,2);
elseif ks ==-1% 原资料阵
    X = X;
end
%% 基本处理
[m,n] = size(X);mnl = min(m,n);
if m<n %无需时间转换
    S = X*X'/(n-1);
    [lambda,L] = jacobi_eig(S,10E7);
    Y = L'*X;
else %需时间转换
    S_Q= X'*X/(m-1);
    [lambda_Q,L_Q] = jacobi_eig(S_Q,10E7);
    for i = 1:mnl
        L(:,i) = X*L_Q(:,i) /sqrt((m-1)*lambda_Q(i));
        lambda(i) = (m-1)*lambda_Q(i)/(n-1);
    end
    Y =L'*X;
end

%% 赋物理意义
if kks ==1 % 选择赋物理意义
    for i = 1:length(lambda)
        L(:,i) = L(:,i)*sqrt(lambda(i));
        Y(i,:) = Y(i,:)/sqrt(lambda(i));
    end
end

%% 计算方差贡献率
vf = zeros(length(lambda),1);cvf = zeros(length(lambda),1);
for i = 1:length(lambda)
    vf(i) = lambda(i)/sum(lambda);
end
for i = 1:length(lambda)
    cvf(i) =sum(vf(1:i));
end

%% 辅助函数
function [a,L] = jacobi_eig(A,maxit,tol)
% a的元素为A的特征值，V的第j列为对应于特征值a(j,1)的特征向量
if nargin == 2
    tol = 1.0e-6;
elseif nargin == 1
    maxit = 1000;
    tol = 1.0e-6;
end
A0 = A;
n = size(A,1);
V = eye(n);
b = A0 - diag(diag(A0));
w = sqrt(sum(sum(b .* b)));
[r,c] = find(abs(b) > w);
k = 1;
while (w > tol) && (k < maxit)
    if (~isempty(r))
        for i = 1:size(r,1)
            if (r(i,1) > c(i,1))
                if (abs(A0(r(i,1),c(i,1))) < tol)
                    u = 0;
                else
                    u = acot((A0(r(i,1),r(i,1)) - A0(c(i,1),c(i,1)))...
                        / (2 * A0(r(i,1),c(i,1)))) / 2;
                end
                V1 = eye(n);
                V1(r(i,1),r(i,1)) = cos(u);
                V1(c(i,1),c(i,1)) = cos(u);
                V1(r(i,1),c(i,1)) = -sin(u);
                V1(c(i,1),r(i,1)) = sin(u);
                A0 = V1' * A0 * V1;
                V = V * V1;
                k = k + 1;
            end
        end
    end
    b = A0 - diag(diag(A0));
    w = w / n;
    [r,c] = find(abs(b) > w);
end
if (k <= maxit)
    disp(['Jacobi方法迭代',num2str(k),'步收敛！']);
else
    disp(['迭代步数大于',num2str(maxit),'收敛太慢了！']);
end
a = diag(A0);

[~,ind] = sort(a,'descend');
a = a(ind);Vs = V(:,ind);
L =Vs;
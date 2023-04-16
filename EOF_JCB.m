function [L,Y,lambda,vf,cvf]= EOF_JCB(X,ks,kks)
%�������˵��
%XΪԭ���Ͼ���[m n]��mΪ��������nΪ��������
%ksΪ���ϱ������ã�ks = -1����ԭʼ��ks = 0���Ļ���ks = 1���ñ�׼������
%kksΪ�Ƿ���������ѡ�kks = 1��ʾ���������壻kks = 0��ʾ������������
%
%�������˵��
%LΪϵ������[m mnl]��YΪ���ɷ־���[mnl n]��lambdaΪ��������
%vfΪ������ʣ�cvfΪ�ۻ��������

%% �ж����ϱ�������
if ks == 0% ��ƽ��
    X = X - mean(X,2);
elseif ks == 1% ��׼������
    X =  (X - mean(X,2))./std(X,0,2);
elseif ks ==-1% ԭ������
    X = X;
end
%% ��������
[m,n] = size(X);mnl = min(m,n);
if m<n %����ʱ��ת��
    S = X*X'/(n-1);
    [lambda,L] = jacobi_eig(S,10E7);
    Y = L'*X;
else %��ʱ��ת��
    S_Q= X'*X/(m-1);
    [lambda_Q,L_Q] = jacobi_eig(S_Q,10E7);
    for i = 1:mnl
        L(:,i) = X*L_Q(:,i) /sqrt((m-1)*lambda_Q(i));
        lambda(i) = (m-1)*lambda_Q(i)/(n-1);
    end
    Y =L'*X;
end

%% ����������
if kks ==1 % ѡ����������
    for i = 1:length(lambda)
        L(:,i) = L(:,i)*sqrt(lambda(i));
        Y(i,:) = Y(i,:)/sqrt(lambda(i));
    end
end

%% ���㷽�����
vf = zeros(length(lambda),1);cvf = zeros(length(lambda),1);
for i = 1:length(lambda)
    vf(i) = lambda(i)/sum(lambda);
end
for i = 1:length(lambda)
    cvf(i) =sum(vf(1:i));
end

%% ��������
function [a,L] = jacobi_eig(A,maxit,tol)
% a��Ԫ��ΪA������ֵ��V�ĵ�j��Ϊ��Ӧ������ֵa(j,1)����������
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
    disp(['Jacobi��������',num2str(k),'��������']);
else
    disp(['������������',num2str(maxit),'����̫���ˣ�']);
end
a = diag(A0);

[~,ind] = sort(a,'descend');
a = a(ind);Vs = V(:,ind);
L =Vs;
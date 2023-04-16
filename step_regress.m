function step_regress(data)
    %即在这里输入的data矩阵必须都是横资料阵，最后一行是预报量
    r = corr(data');
    [y,n]=size(data);%y 代表预报量所处的位置 n 代表样本容量
    L=0;K=[];t=y-1; %t是预报因子数，K用于统计已经引入的因子序号，L代表有效引入的次数
%起步阶段
    while 1
        V=zeros(t,1);%为各待选因子的方差贡献
        for i=1:t
            V(i,1)=r(i,y)^2/r(i,i);
        end
        V(K,1)=0;%去除上回已经引入因子的方差贡献
        %引入检验
        [M,I]=max(V(:,1),[],'omitnan');%找出方差贡献最大的待选因子
        
        %检验是否可以引入
        F=(V(I,1)/1)/((r(y,y)-V(I,1))/(n-(L+1)-1));
        Fa=finv(0.8,1,(n-(L+1)-1));
        %
        
        if F>Fa %如果通过了显著检验
            r=inv_terse(r,I);%引入该因子
            L=L+1;K=[K,I];%有效引入次数+1
        else
            return %直接结束递归
        end
        disp(strcat("起步阶段：",num2str(L)));
        disp(strcat("选入因子为: X",num2str(I)));
        disp(strcat("最大方差贡献为：",num2str(V(I,1))));
        disp(strcat("统计量F:",num2str(F),"   Fa:",num2str(Fa)));
        %标准化的回归系数就是r矩阵的最后一列
        b=zeros(y,1);
        for i=1:y-1
            b(i)=r(i,y)*std(data(y,:))/std(data(i,:));
        end
        
        disp(strcat("各引入因子系数为: ",num2str(b(K,1)'),"  其中各系数排列顺序为: ",num2str(K)));
        
        b(y,1)=mean(data(y,:))-b(K,1)'*mean(data(K,:),2);
        
        disp(strcat("截距b0为：",num2str(b(y,1))));
    
        %复相关系数
        R=sqrt(1-r(y,y));
        disp(strcat("复相关系数为：",num2str(R)))
        
        %剩余平方和
        %标准化回归方程的剩余平方和：ryy
        %原始回归方程的剩余平方和：
        yd=data(y,:)-mean(data(y,:),2);
        Q=yd*yd'*r(y,y);
        disp(strcat("剩余平方和为：",num2str(Q)))
        %剩余方差
        d=Q/(n-L-1);
        disp(strcat("剩余标准差差为：",num2str(sqrt(d))))
        disp("--------------------------------------------------------------------------------------------------")
        if L==3
           break %起步阶段结束，开始剔除
        end
    end
%剔除与引入检验
    while 1
        
%剔除检验
        while 1
            V=zeros(t,1);%为各待选因子的方差贡献
            for i=1:t
                V(i,1)=r(i,y)^2/r(i,i);
            end
            %剔除检验：寻找方差贡献最少的因子
            [M,I]=min(V(K,1),[],'omitnan');
            
            %对改因子进行显著性检验
            F=(V(K(I),1)/1)/(r(y,y)/(n-L-1));
            Fa=finv(0.8,1,n-L-1);
            %
            if F>Fa
                disp(strcat("剔除检验："));
                disp(strcat("选入因子为: X",num2str(K(I))));
                disp(strcat("选入因子通过假设检验，不剔除"));
                disp("--------------------------------------------------------------------------------------------------")
                break %意为无法剔除 跳出剔除检验，剔除检验后直接进入引入检验
            else
                disp(strcat("剔除检验："));
                disp(strcat("选入因子为: X",num2str(K(I))));
                disp(strcat("最大方差贡献为：",num2str(V(I,1))));
                disp(strcat("统计量F:",num2str(F),"   Fa:",num2str(Fa)));
                r=inv_terse(r,K(I));%进行二次求逆，可以剔除该因子的引入对r阵造成的影响
                L=L-1;K=K(~ismember(K,K(I)));
            end
            b=zeros(y,1);
            for i=1:y-1
                b(i)=r(i,y)*std(data(y,:))/std(data(i,:));
            end
            disp(strcat("各引入因子系数为: ",num2str(b(K,1)'),"  其中各系数排列顺序为: ",num2str(K)));
            b(y,1)=mean(data(y,:))-b(K,1)'*mean(data(K,:),2);
            disp(strcat("截距b0为：",num2str(b(y,1))));
            
            %复相关系数
            R=sqrt(1-r(y,y));
            disp(strcat("复相关系数为：",num2str(R)))
    
            
            %剩余平方和
            %标准化回归方程的剩余平方和：ryy
            %原始回归方程的剩余平方和：
            yd=data(y,:)-mean(data(y,:),2);
            Q=yd*yd'*r(y,y);
            disp(strcat("剩余平方和为：",num2str(Q)))
            
            %剩余方差
            d=Q/(n-L-1);
            disp(strcat("剩余标准差差为：",num2str(sqrt(d))))
            disp("--------------------------------------------------------------------------------------------------")
        end        
%引入检验
        while 1
            flag=0;
            V=zeros(t,1);
            for i=1:t
                V(i,1)=r(i,y)^2/r(i,i);
            end
            V(K,1)=0;%去除已选择因子
            %引入检验
            [M,I]=max(V(:,1),[],'omitnan');
    
            F=(V(I,1)/1)/((r(y,y)-V(I,1))/(n-(L+1)-1));
            Fa=finv(0.8,1,n-(L+1)-1);
            if F>Fa
                disp(strcat("引入检验："));
                disp(strcat("选入因子为: X",num2str(I)));
                disp(strcat("最大方差贡献为：",num2str(V(I,1))));
                disp(strcat("统计量F:",num2str(F),"   Fa:",num2str(Fa)));
                r=inv_terse(r,I);
                L=L+1;K=[K,I];
                b=zeros(y,1);
                for i=1:y-1
                    b(i)=r(i,y)*std(data(y,:))/std(data(i,:));
                end
                disp(strcat("各引入因子系数为: ",num2str(b(K,1)'),"  其中各系数排列顺序为: ",num2str(K)));
                b(y,1)=mean(data(y,:))-b(K,1)'*mean(data(K,:),2);
                disp(strcat("截距b0为：",num2str(b(y,1))));
                
                %复相关系数
                R=sqrt(1-r(y,y));
                disp(strcat("复相关系数为：",num2str(R)))
    
                %剩余平方和
                %标准化回归方程的剩余平方和：ryy
                %原始回归方程的剩余平方和：
                yd=data(y,:)-mean(data(y,:),2);
                Q=yd*yd'*r(y,y);
                disp(strcat("剩余平方和为：",num2str(Q)))
                
                %剩余方差
                d=Q/(n-L-1);
                disp(strcat("剩余标准差差为：",num2str(sqrt(d))))
                disp("--------------------------------------------------------------------------------------------------")
                break %引入成功，回到剔除
            else
                flag=1;
                disp(strcat("引入检验："));
                disp(strcat("选入因子为: X",num2str(I)));
                disp(strcat("选入因子没有通过假设检验，不引入"));
                disp("--------------------------------------------------------------------------------------------------")
                break
            end
        end
        if flag==1
           break 
        end
       
    end
%最终结果
    b=zeros(y,1);
    for i=1:y-1
        b(i)=r(i,y)*std(data(y,:))/std(data(i,:));
    end
    b(y,1)=mean(data(y,:))-b(K,1)'*mean(data(K,:),2);
    
    %复相关系数
    R=sqrt(1-r(y,y));
    
    %剩余平方和
    %标准化回归方程的剩余平方和：ryy
    %原始回归方程的剩余平方和：
    yd=data(y,:)-mean(data(y,:),2);
    Q=yd*yd'*r(y,y);
    disp(strcat("剩余平方和为：",num2str(Q)))
    %剩余方差
    d=Q/(n-L-1);
    %F统计量
    F=((1-r(y,y))/L)/(r(y,y)/(n-L-1));
    
    disp("-------------------------------------------最终结果-----------------------------------------------")
    disp(strcat("各引入因子系数为: ",num2str(b(K,1)'),"  其中各系数排列顺序为: ",num2str(K)));
    disp(strcat("统计量F:",num2str(F)));
    disp(strcat("截距b0为：",num2str(b(y,1))));
    disp(strcat("复相关系数为：",num2str(R)));
    disp(strcat("剩余平方和为：",num2str(Q)));
    disp(strcat("剩余标准差为：",num2str(sqrt(d))));
end

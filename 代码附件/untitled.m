%%% 建立图 %%%

%读入数据
data=load('附件1.txt');
n=50; %居民点数量
m=3; %仓库数量
d=data(1:n,1:n); %距离矩阵
demand=data(n+1:end,:); %需求矩阵

%建立完全图
G=graph;
for i=1:n
    for j=1:m
        G=addnode(G,['v',num2str(i),',',num2str(j)],'demand',demand(i,j));
    end
end
for j=1:m
    G=addnode(G,['w',num2str(j)],'demand',inf);
end
G=addnode(G,'s');
G=addnode(G,'t');
for i=1:n
    for j=1:m
        G=addedge(G,'s',['v',num2str(i),',',num2str(j)],d(i,j));
        G=addedge(G,['v',num2str(i),',',num2str(j)],'t',0);
    end
end
for j=1:m
    G=addedge(G,'s',['w',num2str(j)],0);
    G=addedge(G,['w',num2str(j)],'t',0);
end
for i=1:n
    for j=1:m
        for k=1:m
            if j~=k
                G=addedge(G,['v',num2str(i),',',num2str(j)],['v',num2str(i),',',num2str(k)],d(i,j)+d(j,k));
            end
        end
    end
end

%求解MST
T=graphminspantree(G,'s');

%%% 转化为插头问题 %%%

%计算每个插座到每个插头的距离
a=zeros(n,m);
for i=1:n
    for j=1:m
        a(i,j)=min([d(i,1)+d(1,j),d(i,2)+d(2,j),d(i,3)+d(3,j)]);
    end
end

%建立插头问题
n2=n*m; %插座数量
m2=m; %插头数量
d2=zeros(n2,m2); %距离矩阵
for i=1:n
    for j=1:m
        for k=1:m2
            d2((i-1)*m+j,k)=d(i,k);
        end
    end
end
for i=1:n2
    for j=1:m2
        d2(i,j)=a(ceil(i/m),(j-1)/m+1)+d2(i,j);
    end
end

%%% 求解插头问题 %%%

%设置求解参数
params.outputflag=0;
params.method=2;
params.TimeLimit=1800;

%求解线性规划
m=1000; %足够大的常数
vtype=repmat('C',m2*n2,1); %变量类型
ctype=repmat('U',2*n2,1); %约束类型
lb=zeros(m2*n2,1); %下界
ub=m*ones(m2*n2,1); %上界
f=zeros(m2*n2,1); %目标系数
A=zeros(2*n2,m2*n2); %约束系数
b=zeros(2*n2,1); %约束值
for i=1:n2
    for j=1:m2
        f((i-1)*m2+j)=d2(i,j);
        A(i,(i-1)*m2+j)=1;
        A(n2+i,(i-1)*m2+j)=1;
        A(i,(j-1)*n2+1:i+1:(j*n2))=-1;
        A(n2+i,(j-1)*n2+1:i+1:(j*n2))=1;
        b(i)=demand(ceil(i/m),(j-1)/m+1);
        b(n2+i)=0;
    end
end
[x,fval,exitflag,output]=gurobi(f,A,b,[],[],lb,ub,ctype,vtype,params);

%%% 整合结果 %%%

%计算运输路径和完成时间
selected=find(abs(x)>1e-4);
routes={};
for i=1:length(selected)
    [row,col]=ind2sub([n2,m2],selected(i));
    if col<=n
        routes{row}=[routes{row},['v',num2str(col),',',num2str(mod(row-1,m)+1)]];
    else
        routes{row}=[routes{row},['w',num2str(mod(row-1,m)+1)]];
    end
end
T2=zeros(n2,n2);
p=1:length(selected);
T2(sub2ind([n2,n2],mod(selected-1,n2)+1,ceil(selected/n2)))=1;
T2=graph(T2);

%映射回原图
G2=graph;
G2=addnode(G2,'s');
G2=addnode(G2,'t');
for j=1:m
    G2=addnode(G2,['w',num2str(j)],'demand',inf);
end
for i=1:n
    G2=addnode(G2,['v',num2str(i)],'demand',sum(demand(i,:)));
end
for i=1:n
    for j=1:m
        G2=addedge(G2,'s',['v',num2str(i),',',num2str(j)],0);
    end
end
for i=1:n
    for j=1:m
        G2=addedge(G2,['v',num2str(i),',',num2str(j)],'t',0);
    end
end
for j=1:m
    for i=1:n
        G2=addedge(G2,['w',num2str(j)],['v',num2str(i)],d(i,j));
    end
end
for i=1:n
    for j=1:m-1
        G2=addedge(G2,['v',num2str(i),',',num2str(j)],['v',num2str(i),',',num2str(j+1)],0);
    end
    G2=addedge(G2,['v',num2str(i),',',num2str(m)],['v',num2str(i),',',num2str(1)],0);
end
for j=1:m
    for k=1:length(routes)
        if any(strfind(routes{k},['w',num2str(j)]))
            G2=addedge(G2,'s',['v',num2str(k),',',num2str(j)],0);
        end
    end
    Tsp=dijkstra(T2,mod(selected-1,n2)+1,find(strncmp(G.Nodes.Name,'w',1)));
    for k=1:length(routes)
        if any(strfind(routes{k},['w',num2str(j)]))
            G2=addedge(G2,['v',num2str(k),',',num2str(j)],'t',Tsp(k));
        end
    end
end
T2=shortestpath(G2,'s','t','Method','positive');
T2=graph(T2);

%输出结果
fprintf('整体完成时间：%f\n',T2.Edges.Weight(T2.findedge('s','t')));
for j=1:3
    selected=find(strncmp(G.Nodes.Name,['w',num2str(j)],2));
    fprintf('仓库%d出发的车辆数量：%d\n',j,length(selected));
    fprintf('完成时间：%f\n',(sum(T2.Edges.Weight(T2.inedges(selected)))+sum(T2.Nodes.Demand(selected))));
end
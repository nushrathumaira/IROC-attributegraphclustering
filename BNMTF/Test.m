G = load('D:\AttributedGraph\data\syn\new5clusterGraph_network.txt'); 
model=BNMTF(G,0,1,50,5);
U=model.U;
[m,n]=size(U);
newU=[];
K=0;
for i =1 : n
    if(sum(U(:,i))~=0)
        newU=[newU,U(:,i)];
        K=K+1;
    end
end
labelM = zeros(m,K);
T=0;
for i = 1 : m
    for j = 1 : K
        if(newU(i,j)>T)
            labelM(i,j)=1;
        end
    end    
end
dlmwrite('D:\AttributedGraph\data\syn\BNMTFnew5clusters_label.txt',labelM','delimiter', '\t', 'precision', 8);    

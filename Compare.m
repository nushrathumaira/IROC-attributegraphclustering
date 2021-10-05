A = load('D:\research\AttributedGraph\data\syn\2clusterGraph_network.txt');
F = load('D:\research\AttributedGraph\data\syn\2clusterGraph_feature.txt'); 

A = sparse(A);
F=[1,2;1,1;0,0];
numF = size(F, 2);
numV = size(F, 1);
F = F +1;
inxF =[];
for i = 1 : numF
    V = F(:,i);   
    Index = sparse(numV,max(V));
    for i = 1:numV
        Index(i, V(i)) = 1;
    end
    inxF = [inxF, Index];
end

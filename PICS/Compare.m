A = load('D:\AttributedGraph\data\real\0_network.txt');
F = load('D:\AttributedGraph\data\real\0_feature.txt'); 

A = sparse(A);
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
F = inxF;
[k, l, Nx, Ny, Qx, Qy, Dnz, DnzF, cA, c2, CAt, C2t] = cc_search_duo(A, F,'name');
dlmwrite('D:\AttributedGraph\data\syn\result\PICS_0_label.txt',Qx,'delimiter', '\t', 'precision', 8); 
dlmwrite('D:\AttributedGraph\data\syn\result\PICS_0_feature_label.txt',Qy,'delimiter', '\t', 'precision', 8); 
NewF=[];
for i=1:k
    ind = find(Qx==i);
    ind = ind';
    h=[];
    for j=1:l        
       ix = find(Qy == j);
       h = [h, F(ind,ix)];
    end
    NewF=[NewF;h];
end
NewF = full(NewF);

[nx ny] = size(F);
[Px, Py] = cc_perm_cohesive(k, l, Qx, Qy, F);

F2 = F(Px,Py);
figure;
plot_binary_duo(F,1,3);
plot_cluster_grid(nx, ny, Nx, Ny);

figure;
plot_binary_duo(A,1,1);
plot_cluster_grid(nx, ny, Nx, Ny);
			 			 

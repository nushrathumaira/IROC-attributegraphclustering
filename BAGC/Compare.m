adj_coo = load('D:\AttributedGraph\data\syn\5clusterGraph_network.txt');
attr_tab = load('D:\AttributedGraph\data\syn\5clusterGraph_feature.txt'); 

adj_coo = sparse(adj_coo);
    % set cluster number    
    clt_num  = 5; 
    
    % set iteration number for optimization
    iter_num = 10;
    
	% run bagc
    [modularity,entropy,time] = bagc(adj_coo,attr_tab,clt_num,iter_num);
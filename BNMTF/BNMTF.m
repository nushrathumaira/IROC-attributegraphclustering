function model=BNMTF(G,loss_option,sparse_option,lambda,max_rank,U0,B0)
%G: the adjancency matrix of the network
%loss_option =0: square loss; 
%            =1: KL-divergence
%sparse_loss =0: all elements used; 
%            =1: positive elements used
%lambda: the regularization parameter
%max_rank: the maximum value for the rank (optional)
%U0: the initial guess for U (optional)
%B0: the initial guess for B (optional)
%model: struct variable with learned U and B as components
    N = size(G,1);
    
    if (nargin < 5)
        max_rank = ceil(N/2);
    end

    %remove non-connected nodes
    active_nodes = sum(G)~=0;
    active_node_indices = find(active_nodes);
    non_active_nodes = ~active_nodes;
    non_active_node_indices = find(non_active_nodes);

    if sum(non_active_nodes)>0
        G = G(active_node_indices,active_node_indices);
        fprintf('\nNOTE: the following non-connected nodes are removed from the graphs:\n');
        disp(non_active_node_indices);
        if nargin>=6
            U0=U0(active_node_indices,:);
        end
    end
    n=size(G,1);
    
    %set up the initial W,H matrices
    if (nargin < 6)
        U0 = rand(n,max_rank);
    end
    if (nargin < 7)
        B0 = rand(max_rank);
        B0 = (B0+B0')/2;
    end
    
    %run the BNMTF
    max_iter=10;
    
    if loss_option==0
        [U1, B] = BNMTF_sq(G, U0, B0, lambda, sparse_option, max_iter);
    elseif loss_option==1
        [U1, B] = BNMTF_kl(G, U0, B0, lambda, sparse_option, max_iter, 0.1);
    else
        error('Wrong loss option!');
    end
    
    if sum(non_active_nodes)>0
        U=zeros(N,max_rank);
        U(active_node_indices,:)=U1;
    else
        U=U1;
    end
    model.U=U;
    model.B=B;
    clear U B U1 B1 U0 B0 active_nodes active_node_indices non_active_nodes non_active_node_indices;    
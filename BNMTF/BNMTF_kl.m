function [U, B] = BNMTF_kl(G, U0, B0, lambda, sparse_option, max_iter, epsilon)
%sparse_option=0: all elements used; =1: positive elements used
    n=size(U0,1);
    U=U0;
    B=B0;
    clear U0 B0;
    issymetric=issym(G);
    if sparse_option==1
        W=sparse(real(G>0));
        count_after=full(sum(W,2)');
        count_before=full(sum(W,1));
    elseif sparse_option~=0
        error('Wrong sparsity option');
    end
    for iter=1:max_iter
        %%%Update U
        tmp1=(U+epsilon)*B';
        if issymetric==1
            tmp2=tmp1;
        else
            tmp2=(U+epsilon)*B;
        end
        if sparse_option==0
            coeff_A=repmat(sum(tmp1+tmp2,1),n,1)./(U+epsilon);
            coeff_B=epsilon*(coeff_A-n*repmat(sum(B,1),n,1)-n*repmat(sum(B,2)',n,1))+lambda;
        elseif sparse_option==1
            coeff_A=full(W*tmp1+W'*tmp2)./(U+epsilon);
            coeff_B=epsilon*(coeff_A-count_before'*sum(B,1)-count_after'*sum(B,2)')+lambda;            
        end
        clear tmp1 tmp2;
        tmp1=B*U';
        tmp2=U*tmp1;
        tmp3=full(G./(tmp2+sqrt(eps)));
        coeff_C=U.*((tmp3*tmp1')+(tmp3'*U*B));
        clear tmp1 tmp2 tmp3;
        U=update_U(coeff_A,coeff_B,coeff_C);
        clear coeff_A coeff_B coeff_C;
        %%%Update B
        if sparse_option==0
            tmp=sum(U,1);
            B=(B.*(U'*(G./(U*B*U'+sqrt(eps)))*U))./(tmp'*tmp);
            B=full(B);
            clear tmp;
        else
            B=(B.*(U'*(G./(U*B*U'+sqrt(eps)))*U))./(U'*W*U);
            B=full(B);
        end
    end
    clear G W count_after count_before;
    
function U=update_U(coeff_A,coeff_B,coeff_C)
    [n,k]=size(coeff_A);
    U=zeros(n,k);
    for i=1:n
        for j=1:k
            a=coeff_A(i,j);
            b=coeff_B(i,j);
            c=coeff_C(i,j);
            if a>0
                u1=(-b+sqrt(b^2+4*a*c))/(2*a);
                if u1>=0&&u1<=1
                    U(i,j)=u1;
                else
                    U(i,j)=1;
                end
            elseif b>0
                if c<=b
                    U(i,j)=c/b;
                else
                    U(i,j)=1;
                end
            elseif b<0
                U(i,j)=1;
            else
                U(i,j)=1;
            end
        end
    end
    clear coeff_A coeff_B coeff_C;
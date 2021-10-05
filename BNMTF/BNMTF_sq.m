function [U, B] = BNMTF_sq(G, U0, B0, lambda, sparse_option, max_iter)
%sparse_option=0: all elements used; =1: positive elements used
    [n,k]=size(U0);
    U=U0;
    B=B0;
    clear U0 B0;
    issymetric=issym(G);
    tmp1=full(U*B*U'-G);
    tmp2=U*B;
    tmp3=B*U';
    if sparse_option==0
        used_index=[];
    elseif sparse_option==1
        used_index=find(G>0);
    else
        error('Wrong sparsity option');
    end
    
    for iter=1:max_iter
        %%%Update U
        for p=1:n
            for q=1:k
                %%%Find t
                Z1=sparse(n,n);
                Z1(:,p)=tmp2(:,q);
                Z1(p,:)=Z1(p,:)+tmp3(q,:);
                Z1_index=find(Z1>0);
                z2=B(q,q);
                if isempty(used_index)
                    a=z2^2;
                    b=full(2*z2*Z1(p,p));
                    c=full(2*z2*tmp1(p,p)+norm(Z1,'fro')^2);
                    d=full(2*tmp1(Z1_index)'*Z1(Z1_index)+lambda);
                else
                    if G(p,p)>0
                        a=z2^2;
                    else
                        a=0;
                    end
                    intersect_index=intersect(used_index,Z1_index);
                    b=full(2*sqrt(a)*Z1(p,p));
                    c=full(2*sqrt(a)*tmp1(p,p)+norm(Z1(intersect_index),2)^2);
                    d=full(2*tmp1(intersect_index)'*Z1(intersect_index)+lambda);
                end
                [t,flag]=find_optimal_U(a,b,c,d,[-U(p,q),1-U(p,q)]);
                clear Z1_index a b c d;
                %%%Update tmp1, tmp2, tmp3
                if t~=0
                    tmp1=full(tmp1+t*Z1);
                    tmp1(p,p)=tmp1(p,p)+(t^2)*z2;
                    tmp2(p,:)=tmp2(p,:)+t*B(q,:);
                    tmp3(:,p)=tmp3(:,p)+t*B(:,q);
                end
                clear Z1 z2;
                %%%Update U(p,q)
                if flag==-1
                    U(p,q)=U(p,q)+t;
                elseif flag==0
                    if abs(U(p,q)+t)>sqrt(eps)
                        error('wrong left flag');
                    end
                    U(p,q)=0;
                else
                    if abs(U(p,q)+t-1)>sqrt(eps)
                        error('wrong right flag');
                    end
                    U(p,q)=1;
                end
            end
        end
        %%%Update B
        if issymetric==0
            for p=1:k
                for q=1:k
                    %%%Find t
                    %P1 equals tmp1
                    P0=U(:,p)*U(:,q)';
                    if isempty(used_index)
                        coeff=-U(:,p)'*tmp1*U(:,q)/(norm(U(:,p),2)*norm(U(:,q),2))^2;
                    else
                        coeff=-tmp1(used_index)'*P0(used_index)/(norm(P0(used_index),2)^2);
                    end
                    t=max(-B(p,q),coeff);
                    %%%Update tmp1, tmp2, tmp3;
                    if t~=0
                        tmp1=tmp1+t*P0;
                        tmp2(:,q)=tmp2(:,q)+t*U(:,p);
                        tmp3(p,:)=tmp3(p,:)+t*U(:,q)';
                        if coeff>-B(p,q)
                            B(p,q)=B(p,q)+t;
                        else
                            B(p,q)=0;
                        end
                    end
                    clear P0;
                end
            end
        else
            for p=1:k
                %%%Update B(p,p)
                P0=U(:,p)*U(:,p)';
                if isempty(used_index)
                    coeff=-U(:,p)'*tmp1*U(:,p)/(norm(U(:,p),2)^4);
                else
                    coeff=-tmp1(used_index)'*P0(used_index)/(norm(P0(used_index),2)^2);
                end
                t=max(-B(p,p),coeff);
                %%%Update tmp1, tmp2, tmp3;
                if t~=0
                    tmp1=tmp1+t*P0;
                    tmp2(:,p)=tmp2(:,p)+t*U(:,p);
                    tmp3(p,:)=tmp3(p,:)+t*U(:,p)';
                    if coeff>-B(p,p)
                        B(p,p)=B(p,p)+t;
                    else
                        B(p,p)=0;
                    end
                end
                clear P0;
                for q=p+1:k
                    %%%Update B(p,q) and B(q,p)
                    P2=U(:,p)*U(:,q)'+U(:,q)*U(:,p)';
                    if isempty(used_index)
                        coeff=-(U(:,p)'*tmp1*U(:,q)+U(:,q)'*tmp1*U(:,p))/(2*((U(:,p)'*U(:,q))^2+(norm(U(:,p),2)*norm(U(:,q),2))^2));
                    else
                        coeff=-tmp1(used_index)'*P2(used_index)/(norm(P2(used_index),2)^2);
                    end
                    t=max(-B(p,q),coeff);
                    if t~=0
                        tmp1=tmp1+t*P2;
                        tmp2(:,q)=tmp2(:,q)+t*U(:,p);
                        tmp2(:,p)=tmp2(:,p)+t*U(:,q);
                        tmp3(p,:)=tmp3(p,:)+t*U(:,q)';
                        tmp3(q,:)=tmp3(q,:)+t*U(:,p)';
                        if coeff>-B(p,q)
                            B(p,q)=B(p,q)+t;
                            B(q,p)=B(p,q);
                        else
                            B(p,q)=0;
                            B(q,p)=0;
                        end
                    end
                    clear P2;
                end
            end
        end
    end
    clear tmp1 tmp2 tmp3 G;
    
function [t,flag]=find_optimal_U(a,b,c,d,range)
%flag=-1: corresponding element in U is between two endpoints
%    = 0: corresponding element is the left endpoint
%    = 1: corresponding element is the right endpoint
    flag=-1;
    left=range(1);
    right=range(2);
    if a>0
        r=c/(2*a)-3*b^2/(16*a^2);
        s=b^3/(32*a^3)-b*c/(8*a^2)+d/(4*a);
        Delta=r^3/27+s^2/4;
        if Delta>0
            t1=(sqrt(Delta)-s/2)^(1/3)-(sqrt(Delta)+s/2)^(1/3)-b/(4*a);
            if t1>=right
                t=right;
                flag=1;
            elseif t1<=left
                t=left;
                flag=0;
            else
                t=t1;
            end
        elseif Delta==0
            t1=-2*(s/2)^(1/3)-b/(4*a);
            t2=(s/2)^(1/3)-b/(4*a);
            if s==0
                if t1>=right
                    t=right;
                    flag=1;
                elseif t1<=left
                    t=left;
                    flag=0;
                else
                    t=t1;
                end
            else
                if t1>=left&&t1<=right&&t2>=left&&t2<=right
                    f1=a*(t1^4)+b*(t1^3)+c*(t1^2)+d*t1;
                    f2=a*(t2^4)+b*(t2^3)+c*(t2^2)+d*t2;
                    if f1<f2
                        t=t1;
                    elseif f1>f2
                        t=t2;
                    else
                        t=max(t1,t2);
                    end
                elseif t1>=left&&t1<=right
                    t=t1;
                elseif t2>=left&&t2<=right
                    t=t2;
                elseif t1<left
                    t=left;
                    flag=0;
                else
                    t=right;
                    flag=1;
                end
            end
        else
            rho=sqrt(-r^3/27);
            phi=acos(-s/(2*rho));
            t1=2*(rho^(1/3))*cos(phi/3)-b/(4*a);
            t2=2*(rho^(1/3))*cos((phi+2*pi)/3)-b/(4*a);
            t3=2*(rho^(1/3))*cos((phi+4*pi)/3)-b/(4*a);
            t_array=sort([t1 t2 t3]);
            t1=t_array(1);
            t2=t_array(2);
            t3=t_array(3);
            if t1>=left&&t1<=right&&t3>=left&&t3<=right
                f1=a*(t1^4)+b*(t1^3)+c*(t1^2)+d*t1;
                f3=a*(t3^4)+b*(t3^3)+c*(t3^2)+d*t3;
                if f1<f3
                    t=t1;
                elseif f1>f3
                    t=t3;
                else
                    t=max(t1,t3);
                end
            elseif t1>=left&&t1<=right
                t=t1;
            elseif t3>=left&&t3<=right
                t=t3;
            elseif t2>=left&&t2<=right
                f1=a*(left^4)+b*(left^3)+c*(left^2)+d*left;
                f2=a*(right^4)+b*(right^3)+c*(right^2)+d*right;
                if f1<=f2
                    t=left;
                    flag=0;
                else
                    t=right;
                    flag=1;
                end
            elseif t1>right
                t=right;
                flag=1;
            elseif t3<left
                t=left;
                flag=0;
            elseif t2>right
                t=left;
                flag=0;
            else
                t=right;
                flag=1;
            end
        end
    elseif b>0
        if c^2<3*b*d
            t=left;
            flag=0;
        elseif c^2>3*b*d
            t1=(-c+sqrt(c^2-3*b*d))/(3*b);
            t2=(-c-sqrt(c^2-3*b*d))/(3*b);
            if t1<=right&&t1>=left
                t=t1;
            elseif t2<=right&&t2>=left
                f1=b*(left^3)+c*(left^2)+d*left;
                f2=b*(right^3)+c*(right^2)+d*right;
                if f1<=f2
                    t=left;
                    flag=0;
                else
                    t=right;
                    flag=1;
                end
            elseif t1<left||t2>right
                t=left;
                flag=0;
            else
                t=right;
                flag=1;
            end
        else            
            t1=-c/(3*b);
            if t1==right
                t=t1;
                flag=1;
            elseif t1<right&&t1>left
                t=t1;
            else
                t=left;
                flag=0;
            end
        end
    elseif c~=0
        t1=-d/(2*c);
        if c>0
            if t1<=left
                t=left;
                flag=0;
            elseif t1>=right
                t=right;
                flag=1;
            else
                t=t1;
            end
        else
            if abs(left-t1)>=abs(right-t1)
                t=left;
                flag=0;
            else
                t=right;
                flag=-1;
            end
        end
    elseif d~=0
        if d>0
            t=left;
            flag=0;
        else
            t=right;
            flag=1;
        end
    else
        t=0;
    end
    if flag==-1
        if t==left
            flag=0;
        elseif t==right
            flag=1;
        end
    end
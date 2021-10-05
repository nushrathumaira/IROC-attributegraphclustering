function result=issym(A)
    result=isequal(A,A.');
    %Other candidates:
    % result=~(nnz(A-A')>0);
    % result=all(all(A==A.'));
    % result=(max(max(A-A'))<=sqrt(eps));
function res = chol_solve_diag(a,z)
    %a: the cholesky factorization
    %z: is b in the usual notation inv(a)*b
    if not(size(a,1)==size(z,2))
        error('a and z must be square matrices with the same dimension');
    end
    K=size(a,1);
    res=zeros(K,1);
    b=(a'\z);
    x=zeros(K,1);
    for h=K:-1:1
        for i=K:-1:h
            x(i)=1/a(i,i)*(b(i,h)-a(i,i+1:end)*x(i+1:end));
        end
        res(h)=x(h);
    end
end
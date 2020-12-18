
function [x, relres, it] = jacobi(A,b,tol,maxit)
    
    n = size(A,1);
    x0 = zeros(n,1);
    normb = norm(b);
    
 resvec = zeros(maxit,1);
 
res = b-A*x0;
relres = norm(res)/normb;

Dm1 = 1 / diag(A);

it  = 0;
while (relres > tol) & (it < maxit)
    it = it + 1;
    x = x0 + Dm1*res;
    x0 = x;
    res = b-A*x0;
    relres = norm(res)/normb;
    resvec(it) = relres;
end

endfunction

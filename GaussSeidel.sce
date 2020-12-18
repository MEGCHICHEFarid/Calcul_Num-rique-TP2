

function [x, relres, it] = GaussSeidel(A,b,tol,maxit)
    
    n = size(A,1);
    x0 = zeros(n,1);
    normb = norm(b);
    
 resvec = zeros(maxit,1);
 
res = b-A*x0;
relres = norm(res)/normb;

DmE = 1 / tril(A);

quali = zeros(maxit,1);

it  = 0;
while (relres > tol) & (it < maxit)
    it = it + 1;
    x = DmE\res;
    x0 = x + x0;
    res = b-A*x0;
    relres = norm(res)/normb;
    resvec(it) = relres;
end

endfunction

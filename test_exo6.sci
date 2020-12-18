n = 50
T0=5;
T1= 20;
[A,b,xex] = poisson1D(n,T0,T1);
tol = 1e-8; maxit = 10000;

tic;
[xJ, relresJ, resvecJ, itJ] = jacobi(A,b,tol,maxit);
toc

tic;
[xGS, relresJ, resvecGS, itGS,quali]= gaussSeidel(A,b,tol,maxit);
toc
figure;
plot(log10(resvecJ),'k.');
plot(log10(resvecGS),'r.');

lambdaMin = 4*sin((%pi / 2)*(1/(n+1)))*sin((%pi / 2)*(1/(n+1)));
lambdaMin = 4*sin(n*(%pi / 2)*(1/(n+1)))*sin((%pi / 2)*(1/(n+1)));

alpha = 2/(lambdaMin+lambdaMax)

tic
[xR, relresR, resvecR, itR] = richardson_scalar(A,b,alpha,tol,maxit);
toc 

plot(log10(resvecR), 'o');

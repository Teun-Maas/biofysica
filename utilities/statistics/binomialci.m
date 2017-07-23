function [lb,ub] = binomialci(z,N,alpha)

nu1 = 2*z;
nu2 = 2*(N-z+1);
F   = finv(alpha/2,nu1,nu2);
lb  = (nu1.*F)./(nu2 + nu1.*F);


nu1 = 2*(z+1);
nu2 = 2*(N-z);
F   = finv(1-alpha/2,nu1,nu2);
ub = (nu1.*F)./(nu2 + nu1.*F);
function S = mcsim(S0,T,L,r,sig,n)
epsilon = randn(L,n);
dt = T/L;
S = S0*prod(exp((r - sig^2/2)*dt+sig*sqrt(dt)*epsilon),1);

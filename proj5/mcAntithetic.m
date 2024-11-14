function S = mcAntithetic(S0,r,sig,T,L,n)
epsilon = randn(L,n);
dt = T/L;
S = S0*prod(exp((r - sig^2/2)*dt+sig*sqrt(dt)*[epsilon -epsilon]),1);

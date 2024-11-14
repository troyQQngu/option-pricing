function C0 = mcControlVar(S0,K,r,T,sig,L,n)

epsilon = randn(L,n);
dt = T/L;
Sa = S0*prod(exp((r - sig^2/2)*dt+sig*sqrt(dt)*epsilon),1);
Sb = S0*prod(exp((r - sig^2/2)*dt+sig*sqrt(dt)*-epsilon),1);

payoff_a = max([Sa(end,:)-K;zeros(1,n)],[],1);
payoff_b = max([Sb(end,:)-K;zeros(1,n)],[],1);

C0a = exp(-r*T)*payoff_a;
C0b = exp(-r*T)*payoff_b;

c = -mean((C0a-mean(C0a)).*(C0b-mean(C0b)))/mean((C0b-mean(C0b)).^2);

[Call_b,Put_b]= blsprice(S0,K,r,T,sig);

C0 = mean(C0a)-(mean(C0b)-Call_b);






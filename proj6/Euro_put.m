% Problem b: European put option price
function price = Euro_put(S0,K,r,T,sigma,N,M,Smax)
% S0 - stock price at time zero
% K - strike price
% r - risk-free interest rate
% T - time of maturity
% sigma - volatility
% N - resolution of time
% M - resolution of stock price
% Smax - boundary stock price

dt = T/N;
t = 0:dt:T;
dS = Smax/M;
s = 0:dS:Smax;
P = NaN*ones(N+1,M+1);

P(end,:)=max(K-(0:M)*dS,0);

P(:,1) = K;
P(:,end) = 0;

E = (1:M-1)';
a = r/2*E*dt-1/2*sigma^2*E.^2*dt;
b = 1+sigma^2*E.^2*dt+r*dt;
c = -r/2*E*dt-1/2*sigma^2*E.^2*dt;
D = spdiags([[a(2:end);0] b [0;c(1:end-1)]],-1:1,M-1,M-1);

for i = N:-1:1
    y = P(i+1,2:end-1)'+[-a(1)*K; zeros(M-3,1); -c(end)*0];
    P(i,2:end-1) = D\y;
end

price = interp1(s,P(1,:),S0);

end
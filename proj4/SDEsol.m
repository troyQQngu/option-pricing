function X = SDEsol(X0,mu,sigma,T,N)
dt = T/N;
X = zeros(N+1,1);
a = zeros(N,1);
X(1) = X0;
for i=1:N
    a(i) = randn;
    X(i+1) = X(i)*exp((mu-1/2*sigma^2)*dt+sigma*sqrt(dt)*a(i));
end
end
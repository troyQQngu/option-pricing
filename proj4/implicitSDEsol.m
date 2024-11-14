function X = implicitSDEsol(X0,mu,sigma,theta,N,T)

dt = T/N;

X = zeros(N+1,1);
X(1) = X0;
a = zeros(T,1);
for i = 1:N
    a(i) = randn;
    X(i+1) = (X(i)+(1-theta)*dt*mu*X(i)+sqrt(dt)*a(i)*sigma*X(i))/(1-theta*dt*mu);
end

end
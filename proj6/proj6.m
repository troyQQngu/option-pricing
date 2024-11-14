S0 = 100;
K = 100;
r = 0.02;
sigma = 0.3;
T = 5/12;
Smax = 300;

p_BS = BS_price(S0,K,r,sigma,T);

N = 10:10:100;% resolution of time
n = length(N);
M = 10:10:100;% resolution of stock price
m = length(M);
error = zeros(n,m);

for i = 1:n
    for j = 1:m
        p = Euro_put(S0,K,r,T,sigma,N(i),M(j),Smax);
        error(i,j) = p - p_BS;
    end
end


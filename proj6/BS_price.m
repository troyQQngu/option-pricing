function p = BS_price(S0,K,r,sigma,T)

d1 = (log(S0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T));
d2 = d1-sigma*sqrt(T);

p = exp(-r*T)*K*normcdf(-d2)-S0*normcdf(-d1);

end
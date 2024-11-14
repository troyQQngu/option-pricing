a = 9343;
m = 62311;

U = zeros(100000,1);
U(1) = mod(a,m);

for i = 1:length(U)-1
    U(i+1) = mod(a*U(i),m);
end
U = U/m;

u = U(1:99999);
v = U(2:100000);
figure(3)
plot(u,v,LineStyle='none',Marker='*')


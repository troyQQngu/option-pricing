% MATH 512 Proj 4
clear all
close all
%% 4a
% define parameters
T = 20;
mu = -2;
sigma = 0.2;
theta = 0.2;
N = 100;
dt = T/N;

% simulate results
X_sim = zeros(N+1,1);
X_sim(1) = 2;
a = zeros(T,1);
for i = 1:N
    a(i) = randn;
    X_sim(i+1) = (X_sim(i)+(1-theta)*dt*mu*X_sim(i)+sqrt(dt)*a(i)*sigma*X_sim(i))/(1-theta*dt*mu);
end

figure(1)
plot(0:dt:T,X_sim,LineWidth=2.5)
hold on
%% 4b
% compute the analytical solution
X_ana = zeros(N+1,1);
X_ana(1) = 2;
for i=1:N
    X_ana(i+1) = X_ana(i)*exp((mu-1/2*sigma^2)*dt+sigma*sqrt(dt)*a(i));
end

plot(0:dt:T,X_ana,'--',LineWidth=2.5)
title('Implicit Method solution V.S. Analytical Solution',FontSize=15)
xlabel('t')
ylabel('X(t)')
legend('Implicit Method','Analytical Solution')
%% 4c
% to verify that the condition for mean square stability, we run multiple
% trials to see if X(t) converges

figure(2)
mu = -2;
sigma = 0.2;
T = 20;
N = 100;
dt = T/N;
n = 1000; % number of trials
Xend = zeros(n,1);
subplot(3,1,1)
for i =1:n
    X = SDEsol(2,mu,sigma,T,N);
    plot(0:dt:T,X)
    hold on
    Xend(i)=X(end);
end
title('Analytical Solution Simulation (sigma=0.2 mu=-2, mean square stable condition satisfied)')
xlabel('t')
ylabel('X(t)')
disp(mean(Xend.^2))

% mu = -1/2*sigma^2
mu = -0.02;
sigma = 0.2;

subplot(3,1,2)
for i =1:n
    X = SDEsol(2,mu,sigma,T,N);
    plot(0:dt:T,X)
    hold on
    Xend(i)=X(end);
end
title('Analytical Solution Simulation (sigma=0.2 mu=-0.02, mean square stable boundary value)')
xlabel('t')
ylabel('X(t)')
disp(mean(Xend.^2))

% mu>-1/2*sigma^2
mu = 0;
sigma = 0.2;

subplot(3,1,3)
for i =1:n
    X = SDEsol(2,mu,sigma,T,N);
    plot(0:dt:T,X)
    hold on
    Xend(i)=X(end);
end
title('Analytical Solution Simulation (sigma=0.2 mu=-0.01, mean square stable condition not satisfied)')
xlabel('t')
ylabel('X(t)')
disp(mean(Xend.^2))

%% 4d
% to verify mean square stability for implicit method

figure(3)
mu = -2;
sigma = 0.2;
T = 20;
dt = T/N;
Xend = zeros(n,1);
subplot(3,1,1)
for i =1:20
    X = implicitSDEsol(2,mu,sigma,theta,N,T);
    plot(0:dt:T,X)
    hold on
    Xend(i) = X(end);
end
title('Implicit Method Simulation 20 trials (dt=0.02, condition satisfied)')
xlabel('t')
ylabel('X(t)')

% dt = -2(mu+1/2*sigma^2)/((1-2theta)mu^2)
T = 165;
dt = T/N;
subplot(3,1,2)
for i =1:20
    X = implicitSDEsol(2,mu,sigma,theta,N,T);
    plot(0:dt:T,X)
    hold on
end
title('Implicit Method Simulation 20 trials (dt=1.65, condition boundary value)')
xlabel('t')
ylabel('X(t)')

% dt > -2(mu+1/2*sigma^2)/((1-2theta)mu^2)
T = 200;
dt = T/N;
subplot(3,1,3)
for i =1:20
    X = implicitSDEsol(2,mu,sigma,theta,N,T);
    plot(0:dt:T,X)
    hold on
end
title('Implicit Method Simulation 20 trials (dt=2, condition not satisfied)')
xlabel('t')
ylabel('X(t)')

%% 4e
% To verify asymptotic stability
figure(4)
mu = 0.01;
sigma = 0.2;
T = 4000;% let T goes to infinity while keeping the same step size
N = 20000;
dt = T/N;
n=100;
% record the end value of each simulation
Xend = zeros(n,1);
subplot(3,1,1)
for i =1:n
    X = SDEsol(2,mu,sigma,T,N);
    plot(0:dt:T,X)
    hold on
    Xend(i) = X(end);
end
title('Analytical Solution Simulation(sigma=0.2 mu=0.01, asymptotically stable condition satisfied)',FontSize=15)
xlabel('t')
ylabel('X(t)')
ylim([0 100])
disp(mean(Xend))

mu = 0.02;
sigma = 0.2;
Xend = zeros(20,1);
subplot(3,1,2)
for i =1:n
    X = SDEsol(2,mu,sigma,T,N);
    plot(0:dt:T,X)
    hold on
    Xend(i) = X(end);
end
title('Analytical Solution Simulation (sigma=0.2 mu=0.02, asymptotically stable condition boundary value)',FontSize=15)
xlabel('t')
ylabel('X(t)')
ylim([0 100])
disp(max(Xend))

mu = 0.03;
sigma = 0.2;
Xend = zeros(20,1);
subplot(3,1,3)
for i =1:n
    X = SDEsol(2,mu,sigma,T,N);
    plot(0:dt:T,X)
    hold on
    Xend(i) = X(end);
end
title('Analytical Solution Simulation (sigma=0.2 mu=0.03, asymptotically stable condition not satsfied)',FontSize=15)
xlabel('t')
ylabel('X(t)')
ylim([0 100])
disp(max(Xend))




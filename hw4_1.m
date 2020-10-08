clear all
clc
syms y(t)
eqn = diff(y,t) == -0.5*y;%y' + 0.5y = 0
cond = y(0) == 1; % y(0) = 1
ySol(t) = dsolve(eqn,cond)
y0 = 1;
h1_ExplicitEuler = 1.0; %step size
h2_ExplicitEuler = 4.2; %step size
h1_ImplicitEuler = 1.0; %step size
h2_ImplicitEuler = 4.2; %step size

t1 = 0:h1_ExplicitEuler:20-h1_ExplicitEuler;
t2 = 0:h2_ExplicitEuler:20-h2_ExplicitEuler;



% Euler's Method
% Initial conditions and setup
y1_ExplicitEuler = zeros(size(t1));  % allocate the result y
n = numel(y1_ExplicitEuler);  % the number of y values
y1_ExplicitEuler(1) = 1;
% The loop to solve the DE
for t=1:n-1
    y1_ExplicitEuler(t+1) = y1_ExplicitEuler(t) - h1_ExplicitEuler * ySol(t);
end

% Initial conditions and setup
y2_ExplicitEuler = zeros(size(t2));  % allocate the result y
n = numel(y2_ExplicitEuler);  % the number of y values
y2_ExplicitEuler(1) = 1;
% The loop to solve the DE
for t=1:n-1
    y2_ExplicitEuler(t+1) = y2_ExplicitEuler(t) - h1_ExplicitEuler * ySol(t);
end

% Initial conditions and setup
y1_ImplicitEuler = zeros(size(t1));  % allocate the result y
n = numel(y1_ImplicitEuler);  % the number of y values
% The loop to solve the DE
y1_ImplicitEuler(1) = 1;
for t=1:n-1
    f(t+1) = ySol(t+1);
    y1_ImplicitEuler(t+1) = y1_ImplicitEuler(t) - h1_ExplicitEuler * ySol(t);
end

% Initial conditions and setup
y2_ImplicitEuler = zeros(size(t2));  % allocate the result y
n = numel(y2_ImplicitEuler);  % the number of y values
y2_ImplicitEuler(1) = 1;
% The loop to solve the DE
for t=1:n-1
    y2_ImplicitEuler(t+1) = y2_ImplicitEuler(t) - h1_ExplicitEuler * ySol(t);
end


%plot(t,y,t,y1_ExplicitEuler,t,y2_ExplicitEuler,t,y1_ImplicitEuler,t,y2_ImplicitEuler,'-o')
tiledlayout(3,2); % Requires R2019b or later
nexttile
plot(t1,y1_ExplicitEuler,'-o','Color',[1,0,0])%red
title('Explicit Euler method (h = 1.0)')
xlabel('t')
ylabel('y(t)')
nexttile
plot(t2,y2_ExplicitEuler,'-o','Color',[1,1,0])%yellow
title('Explicit Euler method (h = 4.3)')
xlabel('t')
ylabel('y(t)')
nexttile
plot(t1,y1_ImplicitEuler,'-o','Color',[0,0,1])%blue
title('Implicit Euler method (h = 1.0)')
xlabel('t')
ylabel('y(t)')
nexttile
plot(t2,y2_ImplicitEuler,'-o','Color',[0,1,0])%green
title('Implicit Euler method (h = 4.3)')
xlabel('t')
ylabel('y(t)')
nexttile
%ODE45
t = [0 20];
[t,y] = ode45(@(t,y) -0.5*y, t, y0);
plot(t,y)
title('Extract')
xlabel('t')
ylabel('y(t)')
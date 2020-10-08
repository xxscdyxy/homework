clear all
clc
syms y(t)
eqn = diff(y,t) == y;%y' + 0.5y = 0
cond = y(0) == 1; % y(0) = 1
ySol(t) = dsolve(eqn,cond);
y0 = 1;
h1_ExplicitEuler = 0.1; %step size
t1 = 0:h1_ExplicitEuler:20-h1_ExplicitEuler;

% Euler's Method
% Initial conditions and setup
y1_ExplicitEuler = zeros(size(t1));  % allocate the result y
n = numel(y1_ExplicitEuler);  % the number of y values
% The loop to solve the DE
for t=1:n-1
    y1_ExplicitEuler(t+1) = y1_ExplicitEuler(t) + h1_ExplicitEuler * ySol(t);
end
l1 = plot(t1,y1_ExplicitEuler);hold on


%%%e second-order Runge-Kutta method with the step size h = 0.1
f=@(x,y)y; %Write your f(x,y) function, where dy/dx=f(x,y), x(x0)=y0.
x0=0; %example x0=0
y0=1; %example y0=0.5
xn=20;% where we need to find the value of y 
                                            %example x=2
h=0.1; %example h=0.2
t0 = x0:h:xn;
i = 1;
yt0(i) = y0;
while x0<=xn
    k1=h*f(x0,y0); 
    x1=x0+h; 
    k2=h*f(x1,y0+k1);
    y1=y0+(k1+k2)/2;           
    x0=x1;
    y0=y1;  
    yt0(i+1) = y0;
    i = i+1;
end
l2 = plot(t0,yt0);hold on



%%%fourth-order Runge-Kutta method with the step size h = 0.1
h=0.1;                                             % step size
x = 0:h:20;                                         % Calculates upto y(3)
y = zeros(1,length(x)); 
y(1) = 1;                                          % initial condition
F_xy = @(t,r) exp(t);                    % change the function as you desire
for i=1:(length(x)-1)                              % calculation loop
    k_1 = F_xy(x(i),y(i));
    k_2 = F_xy(x(i)+0.5*h,y(i)+0.5*h*k_1);
    k_3 = F_xy((x(i)+0.5*h),(y(i)+0.5*h*k_2));
    k_4 = F_xy((x(i)+h),(y(i)+k_3*h));
    y(i+1) = y(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
end
t11 = 0:h:20;
yt1 = y;
l3 = plot(t11,yt1); hold on


%%%Leapfrog method with the step size h = 0.1
a = 0;
b = 20;
h = 0.1;
n0 = (b-a)/h;
f = @(t) exp(t);
dx0 = (b-a)/n0;
c0 = zeros(size(n0)-1);
value0 = zeros(size(n0)-1);
value0(1)= 1;
a0 = a+dx0;
c0(1) = a0;
for k=2:n
c0(k+1) = a0+k*dx0;
value0(k+1) = value0(k-1) + 2*h*f(c0(k));
end
l4 = plot(c0,value0); hold on

%%%Adams Bashforth method with the step size h = 0.1
[y t] = abm4(@(x,y)y,0,20,1,20/0.1);
l5 = plot(t,y); hold off
title('methods with the step size h = 0.1')
lgd=legend([l1,l2,l3,l4,l5],'Explicit Euler method', ...
    ' second-order Runge-Kutta',...
   'fourth-order Runge-Kutta',' Leapfrog method',...
   ' Adams Bashforth ','NumColumns',1);

xlabel('t')
ylabel('y(t)')

figure
tiledlayout(3,2); % Requires R2019b or later
nexttile
plot(t1,y1_ExplicitEuler)
title('Explicit Euler method')
xlabel('t')
ylabel('y(t)')
nexttile
plot(t0,yt0)%yellow
title('second-order Runge-Kutta')
xlabel('t')
ylabel('y(t)')
nexttile
plot(t11,yt1)%blue
title('fourth-order Runge-Kutta')
xlabel('t')
ylabel('y(t)')
nexttile
plot(c0,value0)%green
title('Leapfrog method')
xlabel('t')
ylabel('y(t)')
nexttile
%ODE45
plot(t,y)
title('Adams Bashforth')
xlabel('t')
ylabel('y(t)')


function [y t] = abm4(f,a,b,ya,n)
h = (b - a) / n;
h24 = h / 24;

y(1,:) = ya;
t(1) = a;

m = min(3,n);

for i = 1 : m % start-up phase, using Runge-Kutta of order 4
    t(i+1) = t(i) + h;
    s(i,:) = f(t(i), y(i,:));
    s2 = f(t(i) + h / 2, y(i,:) + s(i,:) * h /2);
    s3 = f(t(i) + h / 2, y(i,:) + s2 * h /2);
    s4 = f(t(i+1), y(i,:) + s3 * h);
    y(i+1,:) = y(i,:) + (s(i,:) + s2+s2 + s3+s3 + s4) * h / 6;
end

for i = m + 1 : n % main phase
    s(i,:) = f(t(i), y(i,:));
    y(i+1,:) = y(i,:) + (55 * s(i,:) - 59 * s(i-1,:) + 37 * s(i-2,:) - 9 * s(i-3,:)) * h24; % predictor
    t(i+1) = t(i) + h;
    y(i+1,:) = y(i,:) + (9 * f(t(i+1), y(i+1,:)) + 19 * s(i,:) - 5 * s(i-1,:) + s(i-2,:)) * h24; % corrector
end
end



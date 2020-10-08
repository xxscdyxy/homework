clear all
clc
syms y(t)
eqn = diff(y,t) == -y*(1-y);%y' + y(1-y) = 0
cond = y(0) == 1/2; % y(0) = 1
ySol(t) = dsolve(eqn,cond);
a = 0;
b = 1;
[h1,h2,h3,h4,h5,h6] = deal(1,0.1,0.01,0.001,0.0001,0.00001);
[n0,n1,n2,n3,n4,n5,n6] = deal(20,(b-a)/h1,(b-a)/h2,(b-a)/h3,(b-a)/h4,(b-a)/h5,(b-a)/h6);
f = @(t)(1/(exp(t)+1));
value = (f(a)+f(b))/2;
[dx0,dx1,dx2,dx3,dx4,dx5,dx6] = deal((b-a)/n0,(b-a)/n1,(b-a)/n2,(b-a)/n3,(b-a)/n4,(b-a)/n5,(b-a)/n6);
[c0, c1,c2,c3,c4,c5,c6] = deal(zeros(size(n0)-1),zeros(size(n1)-1),zeros(size(n2)-1),zeros(size(n3)-1),zeros(size(n4)-1),zeros(size(n5)-1),zeros(size(n6)-1));
[value0,value1,value2,value3,value4,value5,value6] = deal(zeros(size(n0)-1),zeros(size(n1)-1),zeros(size(n2)-1),zeros(size(n3)-1),zeros(size(n4)-1),zeros(size(n5)-1),zeros(size(n6)-1));
[value0(1),value1(1),value2(1),value3(1),value4(1),value5(1),value6(1)] = deal(1/2,1/2,1/2,1/2,1/2,1/2,1/2);
[a0,a1,a2,a3,a4,a5,a6] = deal(a+dx0,a+dx1,a+dx2,a+dx3,a+dx4,a+dx5,a+dx6);
[c0(1),c1(1),c2(1),c3(1),c4(1),c5(1),c6(1)] = deal(a0,a1,a2,a3,a4,a5,a6);
for k=1:(n0-1)
c0(k+1) = a0+k*dx0;
value0(k+1) = value0(k) + (1/20)*f(c0(k));
end
for k=1:(n1-1)
c1(k+1) = a1+k*dx1;
value1(k+1) = value1(k) + h1*f(c1(k));
end
for k=1:(n2-1)
c2(k+1) = a2+k*dx2;
value2(k+1) = value2(k) + h2*f(c2(k));
end
for k=1:(n3-1)
c3(k+1) = a3+k*dx3;
value3(k+1) = value3(k) + h3*f(c3(k));
end
for k=1:(n4-1)
c4(k+1) = a4+k*dx4;
value4(k+1) = value4(k) + h4*f(c4(k));
end
for k=1:(n5-1)
c5(k+1) = a5+k*dx5;
value5(k+1) = value5(k) + h5*f(c5(k));
end
for k=1:(n6-1)
c6(k+1) = a6+k*dx6;
value6(k+1) = value6(k) + h6*f(c6(k));
end

l0 = plot(c0,value0); hold on
l1 = plot(c1,value1);hold on
l2 = plot(c2,value2);hold on
l3 = plot(c3,value3);hold on
l4 = plot(c4,value4);hold on
l5 = plot(c5,value5);hold on
l6 = plot(c6,value6);hold off
title('linearized trapezoidal method')
lgd=legend([l0,l1,l2,l3,l4,l5,l6],'exact solution', ...
    'h = 1','h = 0.1',...
   'h = 0.01','h = 0.001',...
   'h = 0.0001','h = 0.00001','NumColumns',1);

xlabel('t')
ylabel('y(t)')



[c10, c11,c12,c13,c14,c15,c16] = deal(zeros(size(n0)-1),zeros(size(n1)-1),zeros(size(n2)-1),zeros(size(n3)-1),zeros(size(n4)-1),zeros(size(n5)-1),zeros(size(n6)-1));
[value10,value11,value12,value13,value14,value15,value16] = deal(zeros(size(n0)-1),zeros(size(n1)-1),zeros(size(n2)-1),zeros(size(n3)-1),zeros(size(n4)-1),zeros(size(n5)-1),zeros(size(n6)-1));
[value10(1),value11(1),value12(1),value13(1),value14(1),value15(1),value16(1)] = deal(1/2,1/2,1/2,1/2,1/2,1/2,1/2);
[a0,a1,a2,a3,a4,a5,a6] = deal(a+dx0,a+dx1,a+dx2,a+dx3,a+dx4,a+dx5,a+dx6);
[c10(1),c11(1),c12(1),c13(1),c14(1),c15(1),c16(1)] = deal(a0,a1,a2,a3,a4,a5,a6);
for k=1:(n0-1)
c10(k+1) = a0+k*dx0;
value10(k+1) = value10(k) + (1/20)*(f(c10(k))*(f(c10(k))+1)+f(c10(k+1))*(f(c10(k+1))+1))/2;
end
for k=1:(n1-1)
c11(k+1) = a1+k*dx1;
value11(k+1) = value11(k) + h1*(f(c11(k))*(f(c11(k))+1)+f(c11(k+1))*(f(c11(k+1))+1))/2;
end
for k=1:(n2-1)
c12(k+1) = a2+k*dx2;
value12(k+1) = value12(k) + h2*(f(c12(k))*(f(c12(k))+1)+f(c12(k+1))*(f(c12(k+1))+1))/2;
end
for k=1:(n3-1)
c13(k+1) = a3+k*dx3;
value13(k+1) = value13(k) + h3*(f(c13(k))*(f(c13(k))+1)+f(c13(k+1))*(f(c13(k+1))+1))/2;
end
for k=1:(n4-1)
c14(k+1) = a4+k*dx4;
value14(k+1) = value14(k) + h4*(f(c14(k))*(f(c14(k))+1)+f(c14(k+1))*(f(c14(k+1))+1))/2;
end
for k=1:(n5-1)
c15(k+1) = a5+k*dx5;
value15(k+1) = value15(k) + h5* (f(c15(k))*(f(c15(k))+1)+f(c15(k+1))*(f(c15(k+1))+1))/2;
end
for k=1:(n6-1)
c16(k+1) = a6+k*dx6;
value16(k+1) = value16(k) + h6* (f(c16(k))*(f(c16(k))+1)+f(c16(k+1))*(f(c16(k+1))+1))/2;
end

figure
l10 = plot(c10,value10); hold on
l11 = plot(c11,value11);hold on
l12 = plot(c12,value12);hold on
l13 = plot(c13,value13);hold on
l14 = plot(c14,value14);hold on
l15 = plot(c15,value15);hold on
l16 = plot(c16,value16);hold off
title('direct trapezoidal method')
lgd=legend([l10,l11,l12,l13,l14,l15,l16],'exact solution', ...
    'h = 1','h = 0.1',...
   'h = 0.01','h = 0.001',...
   'h = 0.0001','h = 0.00001','NumColumns',1);

xlabel('t')
ylabel('y(t)')


figure
l2 = plot(c2,value2); hold on
l12 = plot(c12,value12); hold on
title('compare direct and linearized trapezoidal methods')
lgd=legend([l2,l12],'linearized trapezoidal methods', ...
    'direct trapezoidal methods','NumColumns',1);
xlabel('t')
ylabel('y(t)')


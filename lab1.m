% Part 1 
close all
clear 
clc
%Constants
alpha = 0.1;
a=1/2*[1 1 sqrt(2)]';
%Function handle
fm=@(t,m) cross(a,m)+alpha*(cross(a,(cross(a,m))));
%Initial condition
m(:,1)=[1 0 0]';
t(1)=0;

%Step size
h=1e-1;
T=40;
N=T/h;

%RK3 loop
for i=1:N
    %Update time
    t(i+1)=t(i)+h;
    %Update m
    k1=fm(t(i),m(:,i));
    k2=fm(t(i)+h,m(:,i)+h*k1);
    k3=fm(t(i)+h/2,m(:,i)+h*k1/4+h*k2/4);
    m(:,i+1)=m(:,i)+h/6*(k1+k2+4*k3);
end
plot(t,m)
title('Landau–Lifshitz equation')
xlabel('time (t)')
ylabel('m')
grid on
figure

for i=1:length(m)
    unit_m(:,i)=m(:,i)/norm(m(:,i));
end

plot3(unit_m(1,:),unit_m(2,:),unit_m(3,:))
title('Trajectory of normalized vector')
hold on
plot3(a(1),a(2),a(3),'r*')
axis equal
legend('m/|m|','a')

N=[20 40 80 160 320 640];
T=40;
start_val(:,1)=[1 0 0]';
t=0;
for i=1:length(N)
    mat_temp=RK3(fm,start_val,T,N(i));
    mat_final(:,i)=mat_temp(:,end);
end

for i=1:length(N)-1
    error(:,i)=norm(mat_final(:,i)-mat_final(:,i+1));
    h(i) = T/N(i);
end

figure
loglog(h,error)
hold on
loglog(h,h.^2,h,h.^3,h,h.^4)
title('Order of accuracy')
ylabel('Error')
xlabel('stepsize h')
legend('error','h^2','h^3','h^4')
grid on

A=[-alpha*(a(2)^2+a(3)^2) alpha*a(2)*a(1)-a(3) alpha*a(3)*a(1)+a(2);...
    a(3)+alpha*a(1)*a(2) -alpha*(a(3)^2+a(1))^2 alpha*a(3)*a(2)-a(1);...
    alpha*a(1)*a(3)-a(2) alpha*a(2)*a(3)+a(1) alpha*(a(1)^2+a(2))^2];

eigenvals=eig(A);
[x,y] = meshgrid(-4:0.01:4,-4:0.01:4);
z = x+1i*y;
R3 = 1+z+z.^2/2+z.^3/6;
zlevel = abs(R3);
figure
plot(eigenvals(1),'o')
hold on
plot(eigenvals(3),'o')
hold on
contour(x,y,zlevel,[1 1],'k');
grid on
title('RK3 Stability Region')
legend('Stable point','Unstable point')

%finding h_stab for RK3
lambda=real(eigenvals(1));
p=[-lambda^3/6 -lambda^2/2 -lambda -2];
r=roots(p)

% Part 2 
clear
close all 
clc

%Initial conditions
r=[1.2 0]; 
dr=[0 -1];
my=1/82.45;

%Stepsize
T=70;
N=7e4; %number of iterations 
h=T/N;  %1e-3;
t(1)=0;

m(:,1:2)=r;
m(:,3:4)=dr;

%RK3
for i=1:3
    k1=func(t(i)    ,m(i,:));
    k2=func(t(i)+h  ,m(i,:)+h*k1);
    k3=func(t(i)+h/2,m(i,:)+h*k1/4+h*k2/4);
    m(i+1,:)=m(i,:)+h/6*(k1+k2+4*k3);  
    t(i+1)=t(i)+h;
end

for i=4:N
    m(i+1,:)=m(i,:)+h/24*(55*func(t(i),m(i,:))-59*func(t(i-1),m(i-1,:))+37*func(t(i-2),m(i-2,:))-9*func(t(i-3),m(i-3,:)));
    t(i+1)=t(i)+h;
end

figure
plot3(t,m(:,1),m(:,2))
title('Trajectory of the satellite as a function of t')

figure
plot(m(:,1),m(:,2))
hold on
scatter(1-my,0)
hold on
scatter(-my,0)
legend('Satellite','Moon','Earth')
title('Trajectory of the satellite')
axis equal

h=0.001;
h2=h/2;
TOL=0.25;
T=60; 
t1=0:h:T;
t2=0:h2:T;

%Explicit Euler and RK3
Euler_h(:,1:2)=r;
Euler_h(:,3:4)=dr;
Euler_h2(:,1:2)=r;
Euler_h2(:,3:4)=dr;
Rh(:,1:2)=r;
Rh(:,3:4)=dr;
Rh2(:,1:2)=r;
Rh2(:,3:4)=dr;
N=T/h;
N2=T/h2;
for i=1:T/h
    Euler_h(i+1,:)=Euler_h(i,:)+h*func(t1(i),Euler_h(i,:));
    k1=func(t1(i)    ,Rh(i,:));
    k2=func(t1(i)+h  ,Rh(i,:)+h*k1);
    k3=func(t1(i)+h/2,Rh(i,:)+h*k1/4+h*k2/4);
    Rh(i+1,:)=Rh(i,:)+h/6*(k1+k2+4*k3);  
end

for i=1:T/h2
    Euler_h2(i+1,:)=Euler_h2(i,:)+h2*func(t2(i),Euler_h2(i,:));
    
    k1=func(t2(i)     ,Rh2(i,:));
    k2=func(t2(i)+h2  ,Rh2(i,:)+h2*k1);
    k3=func(t2(i)+h2/2,Rh2(i,:)+h2*k1/4+h2*k2/4);
    Rh2(i+1,:)=Rh2(i,:)+h2/6*(k1+k2+4*k3);    
end

r_eulerh2=Euler_h2(:,1:2);
r_eulerh2=r_eulerh2(1:2:end,:);
r_euler=Euler_h(:,1:2);
error=0;
i=1;
while error<=TOL
    error=norm(r_euler(i,:)-r_eulerh2(i,:));
    i=i+1;
end
i
disp(['T_acc for explicit Euler is ' num2str(t(i)) ' seconds.'])
figure
plot(Euler_h(1:i,1),Euler_h(1:i,2))
title('Euler')
r_RK3h2=Rh2(:,1:2);
r_RK3h2=r_RK3h2(1:2:end,:);
r_RK3h=Rh(:,1:2);
error=0;
j=1;
while error<=TOL
    error=norm(r_RK3h(j,:)-r_RK3h2(j,:));
    j=j+1;
end
j
disp(['T_acc for Runge Kutta is ' num2str(t(j)) ' seconds.'])
figure
plot(Rh(1:j,1),Rh(1:j,2))
title('RK3')

%Adams Bashforth 
AB(:,1:2)=r;
AB(:,3:4)=dr;
AB2(:,1:2)=r;
AB2(:,3:4)=dr;
for i=1:3
    k1=func(t1(i)    ,AB(i,:));
    k2=func(t1(i)+h  ,AB(i,:)+h*k1);
    k3=func(t1(i)+h/2,AB(i,:)+h*k1/4+h*k2/4);
    AB(i+1,:)=AB(i,:)+h/6*(k1+k2+4*k3);   
    
    k1=func(t2(i)    ,AB2(i,:));
    k2=func(t2(i)+h2 ,AB2(i,:)+h2*k1);
    k3=func(t2(i)+h2/2,AB2(i,:)+h2*k1/4+h2*k2/4);
    AB2(i+1,:)=AB2(i,:)+h2/6*(k1+k2+4*k3);    
end
for i=4:N
    AB(i+1,:)=AB(i,:)+h/24*(55*func(t1(i),AB(i,:))-59*func(t1(i-1),AB(i-1,:))+37*func(t1(i-2),AB(i-2,:))-9*func(t1(i-3),AB(i-3,:)));
end
for i=4:N2
    AB2(i+1,:)=AB2(i,:)+h2/24*(55*func(t2(i),AB2(i,:))-59*func(t2(i-1),AB2(i-1,:))+37*func(t2(i-2),AB2(i-2,:))-9*func(t2(i-3),AB2(i-3,:)));
end
r_ABh2=AB2(:,1:2);
r_ABh2=r_ABh2(1:2:end,:);
r_ABh=AB(:,1:2);
error=0;
k=1;

while error<=TOL
    error=norm(r_ABh2(k,:)-r_ABh(k,:));
    k=k+1;
end
k
disp(['T_acc for Adams Bashforth is ' num2str(t(k)) ' seconds.'])
figure
plot(AB(1:k,1),AB(1:k,2))
title('Adams Bashforth')

initCond=[r dr];
Tacc=t1(k);
timespan=[0 Tacc];
opts=odeset('RelTol',1e-4);
[t,results]=ode23(@(t,u) (func(t,u'))',timespan,initCond,opts);
figure
plot(results(:,1),results(:,2))
xlabel('x(t)')
ylabel('y(t)')
legend('Satellite')
title('ODE23')
h=diff(t);
disp(['The largest timestep for ODE23 is ' num2str(max(h)) ' seconds.'])
disp(['The smallest timestep for ODE23 is ' num2str(min(h)) ' seconds.'])
t(end)=[];

figure
plot(t,h,'-o')
xlabel('Time (t)')
ylabel('Stepsize (h)')
title('Change of stepsize over time')

%% Part 3

close all
clear
clc

r1=0.04;
r2=1e4;
r3=3e7;
Robertson=@(t,x) [-r1*x(1)+r2*x(2)*x(3); r1*x(1)-r2*x(2)*x(3)-r3*x(2)^2; r3*x(2)^2];
N=[125,250,500,1000,2000];
start=[1;0;0];
T=1;
mat125=RK3(Robertson,start,T,N(1));
mat250=RK3(Robertson,start,T,N(2));
mat500=RK3(Robertson,start,T,N(3));
mat1000=RK3(Robertson,start,T,N(4));
mat2000=RK3(Robertson,start,T,N(5));
t=linspace(0,T,1000);
disp(['The smallest stepsize to obtain a stable solution is ' num2str(T/1000) ' steps.'])
mat1000(:,end)=[];
loglog(t,mat1000(1,:),t,mat1000(2,:),t,mat1000(3,:))
xlabel('Time (t)')   
ylabel('x')
title('Robertsons Problem with smallest stepsize')
legend('x_1','x_2','x_3')
%figure
%plot3(mat1000(1,:),mat1000(2,:),mat1000(3,:))
timespan=[0 1];
initCond=[1 0 0];

opts=odeset('RelTol',1e-3,'AbsTol',1e-3/1000);
[t1,result1]=ode23(Robertson,timespan,initCond,opts);
result1=result1';
disp(['Number of stepsizes for RelTol=10^-3 is ' num2str(length(t1)) '.'  ])

opts=odeset('RelTol',1e-4,'AbsTol',1e-4/1000);
[t2,result2]=ode23(Robertson,timespan,initCond,opts);
result2=result2';
disp(['Number of stepsizes for RelTol=10^-4 is ' num2str(length(t2)) '.'  ])

opts=odeset('RelTol',1e-5,'AbsTol',1e-5/1000);
[t3,result3]=ode23(Robertson,timespan,initCond,opts);
result3=result3';
disp(['Number of stepsizes for RelTol=10^-5 is ' num2str(length(t3)) '.'  ])

opts=odeset('RelTol',1e-11,'AbsTol',1e-11/1000);
[t4,result4]=ode23(Robertson,timespan,initCond,opts);
result4=result4';
disp(['Number of stepsizes for RelTol=10^-11 is ' num2str(length(t4)) '.'  ])


h5=diff(t3);
h11=diff(t4);
t3(end)=[];
t4(end)=[];
figure
plot(t3,h5)
hold on
plot(t4,h11)
title('Change of stepsizes over time')
xlabel('Time (t)')
ylabel('Stepsize')
legend('RelTol=1e-5','RelTol=1e-11')

timespan=[0 1000];

opts=odeset('RelTol',1e-3,'AbsTol',1e-3/1000);
[t5,result5]=ode23s(Robertson,timespan,initCond,opts);

opts=odeset('RelTol',1e-4,'AbsTol',1e-4/1000);
[t6,result6]=ode23s(Robertson,timespan,initCond,opts);

opts=odeset('RelTol',1e-5,'AbsTol',1e-5/1000);
[t7,result7]=ode23s(Robertson,timespan,initCond,opts);
disp(['Number of steps for ODE23s RelTol 1e-5 is ' num2str(length(t7))])

opts=odeset('RelTol',1e-11,'AbsTol',1e-11/1000);
[t8,result8]=ode23s(Robertson,timespan,initCond,opts);
disp(['Number of steps for ODE23s RelTol 1e-11 is ' num2str(length(t8))])

t77=t7(1:end-1);
t88=t8(1:end-1);
figure
plot(t77,diff(t7))
hold on
plot(t88,diff(t8))
legend('RelTol=1e-5','RelTol=1e-11')
title('ODE23s')
figure
plot(t88,diff(t8))
title('RelTol 1e^-11, ODE23s')








































clear
close all
clc

%Constants
global L a b Q0 alpha0 Tout T0 k

L=1;
a=0.2;
b=0.3;
Q0=4000;
alpha0=100;
Tout=20;
T0=50;
k=1;
v=10;

h=[0.05 0.025 0.0125 0.00625];

for i=1:length(h)
    T{i}=A_mat(h(i),v)\f_vek(h(i),v);
    N=L/h(i);
    gridpts{i}=linspace(0,L,N);
    
end

for i = 1:length(T)
    endpts(i)=T{i}(end);
end


plot(gridpts{1},T{1},gridpts{2},T{2},gridpts{3},T{3},gridpts{4},T{4})
xlabel('z')
ylabel('Temperature')
title('Change of temperature')
legend('h=0.05','h=0.025','h=0.0125','h=0.00625')
grid on
for i=1:length(h)-1
    error(i)=abs(endpts(i)-endpts(i+1));
end

figure
loglog(h(1:3),error,h(1:3),h(1:3).^2,h(1:3),h(1:3).^3)
xlabel('Stepsize h')
ylabel('Error')
title('Order of accuracy')
legend('Fel','h^2','h^3')
grid on

v=[1 10 30 100];

for i=1:length(v)
    T{i}=A_mat(h(i),v(i))\f_vek(h(i),v(i));
end

figure
plot(gridpts{1},T{1},gridpts{2},T{2},gridpts{3},T{3},gridpts{4},T{4})
legend('v=1 h=0.05','v=10 h=0.025','v=30 h=0.0125','v=100 h=0.00625')
xlabel('z')
ylabel('Temperature')
title('Change of temperature for different velocities and stepsizes')
grid on

% Part 2

h=0.1; N=49; M=21;

Sn=zeros(N,N);

Sn(1,1:2)=[2 -1];
for i=2:N-1
    Sn(i,[i-1 i i+1])=[-1 2 -1];
end
Sn(N,N-1:N)=[-1 2];

Sm(1,1:2)=[2 -1];
for i=2:M-1
    Sm(i,[i-1 i i+1])=[-1 2 -1];
end
Sm(M,M-1:M)=[-1 2];

%Neumann conditions 
Sm(1,2)=-1+Sm(1,2);
Sm(end,end-1)=-1+Sm(end,end-1);

%Constructing A matrix using Kron
In=eye(N); Im=eye(M);
A=kron(Im,Sn)+kron(Sm,In);
A=1/h^2*A;

%Dirichelet conditions
Qx=zeros(M,N);
Qx(:,1)=40;
Qx(:,end)=400;
Qx=reshape(Qx',[],1);
Qx=1/h^2*Qx;

f=100*ones(length(A),1);
F=f+Qx;
T=A\F;

x=meshgrid(0:h:5);
y=meshgrid(0:h:2);

%Dirichelet conditions
first_row=40*ones(1,M);
last_row=400*ones(1,M);

T = reshape(T,[N,M]);
T=[first_row;T;last_row];
figure
mesh(T)
zlabel('Temperature')
title('Temperature of the block')

disp(['Temperature at (x,y)=(3,1) is ' num2str(T(3.1/h,1.1/h)) '.'])

%% Analytical solution
Tx=@(x) -50*x.^2+322*x+40;

xvec=0:h:5;
yvec=0:h:2;

for i=1:length(xvec)
    for j=1:length(yvec)
        matrice(i,j)=xvec(i);
    end
end

Temp_exact=Tx(matrice);
figure
mesh(Temp_exact)
title('Temperature using the analytical solution')
zlabel('Temperature')

abs_error = abs(Temp_exact-T).^2;
mean_square_error = sum(abs_error(:))/numel(Temp_exact);
disp(['The mean square error between the analytical '...
    'and the numerical is ' num2str(mean_square_error)'.' ])


close all

clear all

clc

format long

%part 2 upggift c

ff=@(x,y) 6000*exp(-5*(x-1)^2-10*(y-1.5)^2);

h=0.1;

for i=1:2

    clear F Qx T

    Sn=[];

    N=(5/h)-1;

    M=(2/h)+1;

    x=0:h:N;

    y=0:h:M;

    Sn=(diag(-ones(N-1,1),1)+diag(2*ones(N,1))+diag(-1*ones(N-1,1),-1));

    Sm=(diag(-ones(M-1,1),1)+diag(2*ones(M,1))+diag(-1*ones(M-1,1),-1));

    Sm(end,end-1)=-2;

    Sm(end,end)=2;

    j=1; k=1;

    for x=(1:N).*h

        k=1;

        for y=(0:M-1).*h

            F(j,k)=ff(x,y); 

            k=k+1;

        end

        j=j+1;

    end

    F=reshape(F,[],1);

    In=eye(N); Im=eye(M);

    A=kron(Im,Sn)+kron(Sm,In);

    A=A*(1/h^2);

    Qx_mat=zeros(M,N);

    Qx_mat(:,1)=40;

    Qx_mat(:,end)=400;

    Qx=Qx_mat/(h^2);

    Qx=reshape(Qx',[],1);

    F1=F+Qx;

    T=A\F1;    

    T=reshape(T,[N,M]);

    T_Value(i)=T(3/h+1,1/h+1);

    h=h/2;

end

 

mesh(T)

% contour(T)

% imagesc  













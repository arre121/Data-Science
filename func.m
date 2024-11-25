function [dzdt] = func(t,m)

%Constants 
my=1/82.45;
N=[0 1;-1 0];

%Initial conditions
r0=[-my 0];
r1=[1-my 0];

z1=m(:,1:2); %position
z2=m(:,3:4); %velocity
dz1dt=z2;
dz2dt=-1*(1-my)*(z1-r0)/norm(z1-r0)^3-my*(z1-r1)/norm(z1-r1)^3+(2*N*z2')'+z1;
dzdt=[dz1dt dz2dt];
end


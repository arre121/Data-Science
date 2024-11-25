%% Assignment 1
set(0,'DefaultFigureVisible','on');
clear all
close all
clc
personalNumber = 980417;
%constants
Lm = 2; 
Rm = 21;
b = 1; 
[J, umax] = lab3robot(personalNumber); 
Kt = 38; 
Km = 0.5;
n = 1/20;
%creating an object 's' 
s= tf('s');
%transferfunction, closed for assignment 1 
G = (n*Kt)/(s*(J*s + b)*(s*Lm + Rm) + Km*Kt*s)
lab3robot(G, personalNumber)


%% Assignment 2
F = 5.835; %F=K 5.835 is ok, K = 114 is unstable, 
%the phase margin rises which rotates the nyquist curve so that's over -1 meaning that it's unstable.   
K = F; 
%create closed loop transfer function
Gc = feedback(F*G,1);   %according to fig 3 F*G, 1 is system 2
%Gc = F*G/(1+F*G)
stepinfo(Gc) 
figure
step(Gc,60)

%% Assignment 3
sys = F*G;  %this is basically the open transfer function
figure
margin(sys) %plot
[Gm,Pm,Wcg,Wcp] = margin(sys);
fb = bandwidth(Gc);
PhaseMargin = Pm 
CrossOverFrequency = Wcp
Bandwidth = fb

%% Assignment 6
newWc = 4*Wcp; %4 times faster specification

%with a new cross over frequency we need to 
%evaluate the new phase margin using evalfr(which gives us the frequency response) and margin
f = newWc
frsp = abs(evalfr(sys,f*1i));   %|G(iw)|

[GmNy,PmNy,WcgNy,WcpNy] = margin(sys/frsp); %PmNy is our new Phase margin

%here we determine the phase compensator
correcFactor = -atand(1/18);
PmNy = -PmNy+Pm-correcFactor; % page 109

%For beta we look at figure 5.13 at page 106. Our new phase margin is
%46.6552 degrees and the corresponding beta is approximately 0.16

%we can verify our beta using the following formula
%beta = (1-sind(PmNy))/(1+sind(PmNy))
%this beta = 0.1579 and it shows that 0.16 is good approximation from
%looking at the graph
beta = 0.17;

Td = 1/(newWc*sqrt(beta));
Ti = 18/newWc; %We chose to use 10 in the numerator because this gives a phase 
%reduction of 5.7 degrees which we used earlier. See page 109. However, we increased
%it to 13 because that gave an overshoot less than 5%. 

%We calculate the error coefficient in order to find gamma. 
e1 = (b*Rm+Kt*Km)/(n*Kt); %from page 97
%gamma = 1/(e1*newWc); %from page 110

%The ramp error should be 0.05 or less. Knowing this we can use formula e1
%on page 97 to find gamma. We let s go to zero for the F lead/lag and get the gamma.
%F * G , from page 111, G from assignment 1. When s goes to zero we get
%this 0.05*K2*Kt*n/(b*Rm+Kt*Km = gamma. 
gamma =  0.0432;

%To choose an improved K we use the formula from page 110: 
%K_new * 1/sqrt(beta)*G(iw_cd) = 1*K_old ==> K_new = sqrt(beta)*K_old/G(iw_cd)
%Here we implement that formula for finding new K
K2 = sqrt(beta)*K/(frsp); % frsp is the frequency response at the frequency newWcd

%Now we can create Flead and Flag
Flead = K2*(Td*s+1)/(beta*Td*s+1);

Flag = (Ti*s+1)/(Ti*s+gamma);

%Now we can create the new open transfer function
Fleadlag = Flead*Flag; 
G_new = G*Fleadlag;

%closed transfer function:
G2 = feedback(G_new,1);
stepinfo(G2)

%% Assignment 8
%sid 130

%The formula for sensitivity functions is found on page 61. 
S=1/(1+F*G);
S_leadlag = 1/(1+Fleadlag*G);
figure
bodemag(S,S_leadlag)
legend('Without lead lag','With lead lag')

%% Assignment 9

deltaG1 = (s+10)/40;
deltaG2 = (s+10)/(4*(s+0.01));
T = 1 - S_leadlag;
figure
bode(1/T,deltaG1,deltaG2)
legend('1/T','\DeltaG1','\DeltaG2')
%Reasoning: According to the robust criteria on page 125, for all w the
%delta G should be less than 1/T(iw). This is evident in the plot. 

%% Assignment 10

A = [0 n 0;0 -b/J Kt/J;0 -Km/Lm -Rm/Lm];
B = [0; 0; 1/Lm];
C = [1 0 0];

%To verify if the system is controllable check page 173. To check if the
%system is observable check page 174. 

%If the system is controllable then det(S) should not equal 0.
%If the system is observable then det(O) should not equal 0.

S_control = [B A*B (A^2)*B];
O = [C;C*A;C*(A^2)];

detControl = det(S_control)
detObeserve = det(O)

%Assignment 11

%Both of the determinants are separated from zero, therefore we can say
%that the system is both controllable and observable.


%% Assignment 12 

%The eigenvalues of the matrix A are not all strictly negative. Therefore 
%we must use state feedback to make the system stabile. The current
%eigenvalues of A are 0, -0.5587 and -10.2270. 
eigen = eig(A)
%We can use the function place for pole placement design. 
poles = [-3 -3+1.15*1i -3-1.15*1i]; %these are the desired poles and
%complex numbers need to be conjugates

L = place(A,B,poles);

%forming new A matrix to check for new eigenvalues
A_closed = A - B*L;
eigenA_closed = eig(A_closed)

%creating new closed loop system
D = 0;
sys_closed = ss(A_closed,B,C,D);

%Here we check if the step response of the new system. The system only
%rises to appoximately to 0.0475 rather than 1. 
figure
step(sys_closed)

%To get a 1 we have to scale the input to get 1. 
Kscale = dcgain(sys_closed); %steady state output
K_new = 1/Kscale; %factor to get step response 1
L0 = K_new;
%Now we create a new closed system with the improved step response
sys_closed_scaled = ss(A_closed,B*K_new,C,D);
figure
%Here we check the step response again and see that the gain is 1. 
step(sys_closed_scaled)

lab3robot(G, F, Fleadlag, A, B, C, L, L0, personalNumber);      
















clc; close all; clear all

global Ai Bi Br C Ci dcg t
%% Insert filter values
Rspk = 4;               %Speaker Resistance
Cf = 1.98*1e-6          %Filter Capacitance
Lind = 1.7684*1e-6      %Filter Inductance
Rind = 15*1e-3;         %Shunt resistance
gain = 25;              %Amplifier Gain

Rin = 1500;             %Input Resistor of Aim

%% Tunning Parameters

Ki = 20000; %Integrator weight
gu = pi; %Damping weight

%% Dynamic Model

% Converting to single-ended parameters
Cf = Cf*2; Rspk = Rspk/2;

% State Space Model
A = [-Rind/Lind,    -1/Lind; 
        1/Cf, -1/(Cf*Rspk)];
B = [gain/Lind; 0];
C = [0 1];

G = ss(A, B, [C;zeros(2)], 0);
dcg = dcgain(G); dcg = dcg(1);

%Plotting Bode plot and Step Response
figure(1)
bode(G)
legend('Plant')
grid on

figure(2)
step(G)
legend('Plant')
grid on

%% Integrator Augmentation

Ai = [A,zeros(2,1);-C/gain,0];
Bi = [B;0];
Ci = [C,0];
Br = [0;0;1];

%% LQR - No integrator

R2 = gu;
R1 = C'*dcg^-2*C*2;
N = -dcg^-1*C';

K1 = lqr(A,B,R1,R2,N);
K = [K1,-Ki];
sys1 = ss(Ai-Bi*K,Br+Bi,[Ci;-K;0,0,-K(3)],[0;1;0])

%Plotting Bode plot and Step Response
% figure(1); hold on
% bode(sys1)
% legend('Plant','LQR')
% hold off

figure(2); hold on
step(sys1)
legend('Plant','LQR')
hold off

%% Controller Components

PI.C = 1.5e-9;
PI.R1 = 1/(PI.C*Ki);
PI

tau_i = PI.C*PI.R1

Kim(1) = K(1)*3.33*2; %Remember to multiply by 2 when implementing for BTL
Kim(2) = K(2)*gain;
LQR.Rk = abs(Rin./Kim(1:2))


return
%%

SIMTIME = 0.00005;
Ts = SIMTIME*1e-4;
t = 0:Ts:SIMTIME;

x = particleswarm(@Control_Designer, 3,[0.01;4*1e5],[5; 1e6]);
x(3) = 1;

R2 = x(1);
R1 = C'*dcg^-2*C*2;
N = -dcg^-1*C';

K1 = lqr(A,B,R1,R2,N);
K = [K1,-x(2)];

sys = ss(Ai-Bi*K,Br+Bi*x(3),[Ci;-K;0,0,-K(3)],[0;x(3);0]);

figure(3)
step(G,t); hold on
step(sys,t); hold off
grid on


%%
function [J] = Control_Designer(Q)

global Ai Bi Br C Ci dcg t

R2 = Q(1);
R1 = C'*dcg^-2*C*2;
N = -dcg^-1*C';

K = lqr(Ai(1:2,1:2),Bi(1:2),R1,R2,N);
K = [K,-Q(2)];

% Forming Closed-Loop System
sys = ss(Ai-Bi*K,Br+Bi,Ci,0);

% Calculating Step Response
y = step(sys,t)/dcgain(sys);

% Calculating the cost
J = sum((1-y).^2);
%J = sum(abs(1-y));
%J = t*abs(1-y);
end
clear global, clear variables, close all, clc

bodeopts = bodeoptions;
bodeopts.FreqUnits = 'Hz';

%% Load Data
load('data_UPS')

%% Closed-loop Control Design

% Defining the desired behavior
tso = 3.5e-3; % Measured from Open-Loop experiment
r12 = exp(-4*Ts/tso);
xp = 0.05;

if r12 > 0.97 
    pd1 = exp(-4*Ts/(tso*(1-xp)))*exp(j*0.1); pd2 = conj(pd1);
else
    pd1 = exp(-4*Ts/(tso*(1-xp))); pd2 = pd1^4;
end

fc1 = 50;
wc1 = 2*pi*fc1*Ta; 
pc1 = exp(-2*pi/5); % Cbar PR+Av pole

[Tdez,kt,z1,w] = td_2nd_order(pd1,pd2,Ta,fo);
Sdez = 1 - Tdez;
figure
bode(Tdez,bodeopts)
title('Tdez')
figure
bode(Sdez,bodeopts)
title('Sdez')


%% VRFT Single-loop

Cbar1 = [tf([1],[1],Ta); tf([0 1 0],[1 -2*cos(wc1) 1],Ta); tf([0 0 1],[1 -2*cos(wc1) 1],Ta); tf([0 1 0],[0 1 -pc1],Ta)];
%Cbar1 = [tf([1],[1],Ta); tf([0 1 0],[1 -2*cos(wc1) 1],Ta); tf([0 0 1],[1 -2*cos(wc1) 1],Ta)];

rho1 = vrft_standard(uk,vok1,vok2,tf(Tdez),Cbar1)
Cz1 = minreal(rho1'*Cbar1,1e-3); zpk(Cz1)

%% Remaining control computations

nCz1 = Cz1.num{1};
dCz1 = Cz1.den{1};

zpk(Cz1)

%% Save
save('Controls_designed')
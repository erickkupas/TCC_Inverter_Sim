clear global, clear variables, close all, clc

bodeopts = bodeoptions;
bodeopts.FreqUnits = 'Hz';

fa = 20e3;
Ta = 1/fa;
fo = 50;

% Defining the desired behavior
tso = 4e-3; % Measured from Open-Loop experiment
r12 = exp(-4*Ta/tso);
xp = 0.2;

if r12 > 0.97 
    pd1 = exp(-4*Ta/(tso*(1-xp)))*exp(j*0.1); pd2 = conj(pd1);
else
    pd1 = exp(-4*Ta/(tso*(1-xp))); pd2 = pd1^4;
end

fc1 = 50;
wc1 = 2*pi*fc1*Ta; 
pc1 = exp(-2*pi/5); % Cbar PR+Av pole

[Tdez,kt,z1,w] = td_2nd_order(pd1,pd2,Ta,fo);

sim_time = 1*(1/fo);
N = (sim_time/Ta);
tk = 0:Ta:(N)*Ta;
uk = square(2*pi*fo*tk)';

yk = lsim(Tdez,uk);
figure(); 
plot(tk,yk); grid on;
figure()
e = uk - yk; plot(tk,e); grid on;

N = 100;
Ms = dd_norms(uk,e,Ta,N,'inf') 

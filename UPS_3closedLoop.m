clear global, clear variables, close all, clc

bodeopts = bodeoptions;
bodeopts.FreqUnits = 'Hz';

%%
load('Controls_designed')

%% UPS Parameters
Vin = 400;      %Input Voltage
Vopk = 127*sqrt(2);     %Output voltage (desired) Vorms = 220;
Po = 600; %    Output power -> 600 W

fo = 50;        %Output frequency -> 60 Hz or 50 Hz
fs = 20e3;      %Switching frequency
fa = fs;        %Sampling frequency = Switching frequency
fc = 1850;      %Cut-off frequency -> 30*fo < fc < fs/10 <=> 1800 < fo < 2000
Ro = Vopk^2/(2*Po); %Nominal Load
zeta = 0.4;    %Damping factor
Cf = 1/(4*pi*zeta*fc*Ro); %Capacitor
Lf = 1/(2*pi*fc)^2/Cf; %Inductor

rL = 15e-3;     %Inductor Losses
rC = 1e-6;      %Capacitor Losses
Rdson = 10e-3;  %Transistor Losses
Vsdf = 0.7;     %Transistor Diode
Cin = 3e-3;  %Input Capacitor

%% Parameters for simulation type
enRect = 0;     % 1 to use Rectifier Input, 0 to use pure CC Vin + dVin
dVin = 100;       % Steps for Input voltage (only for enRect = 0)
LoadType = 0;   % 0 for Linear, 1 for NonLinear

%tload_up = 0.04;   % time for step up at Load
%tload_down = 0.07; % time for step down at Load
tload_up = 0.025;
tload_down = 0.035;

t_dvin_up = 1;  % time for step up at Vin 
t_dvin_down = 1;% time for step down at Vin

if (tload_up > tload_down)
    error('Load disturbance: Must use tload_up < tload_down')
end
if (t_dvin_up > t_dvin_down)
    error('Vin disturbance: Must use t_dvin_down < t_dvin_down')
end

sel_control = 1; %0 for Openloop; 1 Single-loop, 2 CascadeI, 3 CascadeII

%% Remaining computations
V3fin = Vin/(sqrt(2*3));    %Equivalent voltage for 3f Input Voltage
Vorms = Vopk/sqrt(2);       %RMS output voltage
RoLin = Vorms^2/Po;         %Total Linear Load
RoLin1 = Vorms^2/(0.2*Po);  %Linear Load 1 (always on)
RoLin2 = Vorms^2/(0.8*Po);  %Linear Load 2 (for steps)
RoNL = Vopk^2/Po;           %Total NonLinear Load
CLoad = 1/(0.03*RoNL*2*fo); %Capacitor for Load (with desired ripple)
RoNL1 = RoNL;               %NonLinear Load 1 (always on)
RoNL2 = RoNL;               %NonLinear Load 2 (always on)
ma = Vopk/Vin;              %Modulation Index for Half-Bridge

%Simulation Time Parameters
Ta = 1/fa;
Ts = 1/fs;
To = 1/fo;
t_step = Ts/100;
%sim_time = 8*To;
sim_time = 30*To;
N = (sim_time/Ta);
tk = 0:Ta:(N)*Ta;

rk = Vopk*sin(2*pi*fo*tk)';
%rk = (Vopk*square(2*pi*fo*tk))';
%ss_time = 5e-3;
%N_ss = ss_time/Ta;
%B = 1/(N_ss/3); % Como definir?
%rk = Vopk*idinput([N+1],'PRBS',[0 B],[-1 1]); % U = idinput(N,TYPE,BAND,LEVELS)
rk = [tk' rk];

fprintf('Simulation parameters OK \n')

%% Salva dados txt
listVars2 = who;
for iVar = 1:length(listVars2)
  S_vars2.(listVars2{iVar}) = eval(listVars2{iVar});
end
save_param_file_psim(S_vars2,'UPS_parameters.txt')

%% Remaining control computations

nCz1 = Cz1.num{1};
dCz1 = Cz1.den{1};

nCez2 = Cez2.num{1};
dCez2 = Cez2.den{1};

nCiz2 = Ciz2.num{1};
dCiz2 = Ciz2.den{1};

nCez3 = Cez3.num{1};
dCez3 = Cez3.den{1};

nCiz3 = Ciz3.num{1};
dCiz3 = Ciz3.den{1};

%% Run Simulink
fprintf('Running simulink with "Control %d" - ',sel_control)
switch sel_control
    case 0
        fprintf('Open-loop\n')
    case 1
        fprintf('Single-loop Control\n')
    case 2
        fprintf('Cascade I Control\n')
    case 3
        fprintf('Cascade II Control\n')
    otherwise
        error('Select control "sel_control" must be an integer between 0 and 3')
end

sim('UPS_sim')

%% Norm
vod = lsim(Tdez,rk(:,2),rk(:,1));
e1 = vod-vok;
Jmr1 = 1/N*e1'*e1


%% Plots

subplot(2,1,1)
plot(rk(:,1),vok,'-k','linewidth',1); grid on; hold on;
plot(rk(:,1),rk(:,2),'--b','linewidth',1);
ylim([-300,300]);
%xlim([0.053,0.115]);
xlabel('Tempo (s)')
legend('v_o(k)','r(k)','location','best')
ylabel('Tensão (V)')
set(gca,'fontname','times') 
set(gca,'fontsize',10)

subplot(2,1,2)
plot(rk(:,1),uk,'-k','linewidth',1); grid on;
xlim([0,0.06]);
ylim([-0.7,0.7]);
xlabel('Tempo (s)')
legend('u(k)','location','best')
ylabel('Ação de controle - u(k)')
set(gca,'fontname','times') 
set(gca,'fontsize',10)

data = table(uk,vok,rk);

writetable(data,'PRAv4.txt');

evk = rk(:,2)-vok;

figure
stairs(evk)

figure
t = 0.1e-3:1/fs:0.2;
vd = zeros(1,length(t));
i = 1;

while t(1,i) <= 1e-3 
    vd(1,i) = 100;
    i = i+1;
end

m1 = (35-100)/(log(4e-3/1e-3));
b1 = 100 - m1*log(1e-3);

while t(1,i) <= 4e-3 
    vd(1,i) = m1*log(t(1,i)) + b1;
    i = i+1;
end


m2 = (20-35)/(log(10e-3/4e-3));
b2 = 35 - m2*log(4e-3);

while t(1,i) <= 10e-3 
    vd(1,i) = m2*log(t(1,i)) + b2;
    i = i+1;
end

while t(1,i) <= 20e-3 
    vd(1,i) = 20;
    i = i+1;
end

while t(1,i) <= 40e-3 
    vd(1,i) = 14;
    i = i+1;
end

while t(1,i) <= 60e-3 
    vd(1,i) = 12;
    i = i+1;
end

while t(1,i) <= 100e-3 
    vd(1,i) = 11;
    i = i+1;
end

while t(1,i) < 200e-3 
    vd(1,i) = 10;
    i = i+1;
end

vd(1,end) = 10;

semilogx(t,vd,'-k','linewidth',1.2); grid on; 
hold on;
semilogx(t,-vd,'-k','linewidth',1.2);
ylim([-110 110]);
xlim([1e-4 0.2]);
ylabel('Desvio relativo da tensão de saída [%]')
xlabel('Duração do transiente [s]')
set(gca,'fontname','times') 
set(gca,'fontsize',10)
text(4e-3,40,' 35% ','fontname','times','FontSize',10)
text(10e-3,27,' 20% ','fontname','times','FontSize',10)
text(21e-3,22,' 14% ','fontname','times','FontSize',10)
text(40e-3,20,' 12% ','fontname','times','FontSize',10)
text(70e-3,18,' 11% ','fontname','times','FontSize',10)
text(118e-3,17,' 10% ','fontname','times','FontSize',10)
text(4e-3,-40,' -35% ','fontname','times','FontSize',10)
text(10e-3,-27,' -20% ','fontname','times','FontSize',10)
text(21e-3,-22,' -14% ','fontname','times','FontSize',10)
text(40e-3,-20,' -12% ','fontname','times','FontSize',10)
text(65e-3,-18,'  -11% ','fontname','times','FontSize',10)
text(118e-3,-17,' -10% ','fontname','times','FontSize',10)
t1 = find(abs(rk(:,1) - 0.0351) < 1e-6);
t2 = find(abs(rk(:,1) - 0.235) < 1e-6);

tek = rk;
tek(:,1) = tek(:,1) - 0.035;

% for aux=1:length(vok(:,1)) 
%     if abs(tek(aux,1)) < 4.75e-3
%         ek(aux) = (vok(aux,1)-rk(aux,2)).*100./rk(aux,2);
%     else
%        ek(aux) = 0; 
%     end
% end

vok(:,1) = vok(:,1) + 1000;
rk(:,2) = rk(:,2) + 1000;

ek = (vok(:,1)-rk(:,2)).*100./rk(:,2);
ek = ek*(-Vopk+1000)/-Vopk;

plot(tek(t1:t2,1),ek(t1:t2),'linewidth',1.2);
%figure
%plot(rk(:,1),ek); grid on;
%xlim([0.0351 0.04])


%% Hinf Norm

clc;
r = rk(:,2);
e = r - vok;
Ta = 1/fa;
N = 100;
[Norm,varargout] = dd_norms(r,e,Ta,N,'inf');
Norm
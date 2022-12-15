clear global, clear variables, close all, clc

%% UPS Parameters

Vin = 400;      %Input Voltage
Vopk = 127*sqrt(2);     %Output voltage (desired);
Po = 600; %    Output power 

fo = 50;        %Output frequency 
fs = 20e3;      %Switching frequency
fa = fs;        %Sampling frequency = Switching frequency
fc = 1850;      % Cut-off frequency -> 30*fo < fc < fs/10 <=> 1800 < fo < 2000
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
dVin = 0;       % Steps for Input voltage (only for enRect = 0)
LoadType = 0;   % 0 for Linear, 1 for NonLinear

tload_up = 1;   % time for step up at load
tload_down = 1; % time for step down at load

t_dvin_up = 1;  % time for step up at Vin
t_dvin_down = 1;% time for step down at Vin

if (tload_up > tload_down)
    error('Load disturbance: Must use tload_up < tload_down')
end
if (t_dvin_up > t_dvin_down)
    error('Vin disturbance: Must use t_dvin_down < t_dvin_down')
end

sel_control = 0; %0 for Openloop; 1 Single-loop, 2 CascadeI, 3 CascadeII

%% Remaining computations

V3fin = Vin/(sqrt(2*3));    %Equivalent voltage for 3f Input Voltage
Vorms = Vopk/sqrt(2);       %RMS output voltage
RoLin = Vorms^2/Po;         %Total Linear Load
RoLin1 = Vorms^2/(0.2*Po);  %Linear Load 1 (always on)
RoLin2 = Vorms^2/(0.8*Po);  %Linear Load 2 (for steps)
RoNL = Vopk^2/Po;           %Total NonLinear Load
CLoad = 1/(0.03*RoNL*2*fo); %Capacitor for Load (with desired ripple)
RoNL1 = RoNL/0.25;          %NonLinear Load 1 (always on)
RoNL2 = RoNL/0.75;          %NonLinear Load 2 (always on)
ma = (Vopk)/Vin;            %Modulation Index for Full-Bridge

%Simulation Time Parameters

%Settling time - Square Wave

Ta = 1/fa;                  %Sampling Time
Ts = 1/fs;                  %Switching Time 
To = 1/fo;                  %Output Period
t_step = Ts/1000;           %Step Time
sim_time = 3*To;            %Simulation Time
Nr = 1448;                  %Samples -> sim_time
N_zeros = 1000;
No = 100;
Nf = N_zeros+No+Nr;
tk = 0:Ta:(Nf)*Ta;
N = No+Nr;
sim_time = Nf*Ta;
trk = 0:Ta:(Nr)*Ta;
%rk_st = Vopk*square(2*pi*fo*trk)';
%rk_st = [zeros(N_zeros+No,1); rk_st];
%rk_st = [tk' rk_st];

% Data collect - PRBS (Pseudorandom binary sequence)

ss_time = 3.5e-3;
N_ss = ss_time/Ta;
B = 1/(N_ss/3); % Como definir?
rk = Vopk*idinput([Nr+1],'PRBS',[0 B],[-1 1]); % U = idinput(N,TYPE,BAND,LEVELS)
length(rk)
%rk = 30*square(2*pi*fo*trk)';
L = length(rk);
fft_rk = fft(rk);
fft_rk = fft_rk(1:L/2);
m_rk = abs(fft_rk)/(L/2);
f = fs*(0:L/2-1)/L;
figure();
plot(f,m_rk,'linewidth',1); grid on;
xlabel('Frequência (Hz)');
ylabel('Amplitude'); 
set(gca,'fontname','times') 
set(gca,'fontsize',10)

rk = [zeros(N_zeros+No,1); rk];
rk = [tk' rk];

fprintf('Simulation parameters OK \n')

%% Salva dados txt
listVars2 = who;
for iVar = 1:length(listVars2)
  S_vars2.(listVars2{iVar}) = eval(listVars2{iVar});
end
save_param_file_psim(S_vars2,'UPS_parameters.txt')

%% Closed-loop Control Design

Cz1 = tf(0,1,Ta);

Cez2 = tf(0,1,Ta);
Ciz2 = tf(0,1,Ta);

Cez3 = tf(0,1,Ta);
Ciz3 = tf(0,1,Ta);

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

vok = vok(N_zeros+1:end-1);
iLk = iLk(N_zeros+1:end-1);
vink = vink(N_zeros+1:end-1);
uk = uk(N_zeros+1:end-1);
vabk = vabk(N_zeros+1:end-1);
iok = iok(N_zeros+1:end-1);

N = length(vok);

%% Noise

SNR = 18^2; %SNRi_dB = 10*log10(SNRi)
iLRMS = sqrt(iLk'*iLk/N); 
voRMS = sqrt(vok'*vok/N);
var_i = iLRMS^2/SNR; var_v = voRMS^2/SNR;

vi1 = sqrt(var_i)*randn(N,1); vi2 = sqrt(var_i)*randn(N,1);
ve1 = sqrt(var_v)*randn(N,1); ve2 = sqrt(var_v)*randn(N,1);

iLk1 = iLk + vi1;
iLk2 = iLk + vi2;
vok1 = vok + ve1;
vok2 = vok + ve2;

%%

save('data_UPS')

%% Plots

figure()
subplot(2,1,1)
%vok = vok(N_zeros:end-1);
rk = rk(N_zeros+1:end-1,1:end);
plot(rk(:,1)*1000,vok+4.26,'linewidth',1); grid on; hold on;
plot(rk(:,1)*1000,rk(:,2),'-k','linewidth',1);
ylim([-3*Vopk,3*Vopk]);
xlim([0.053*1000,0.115*1000]);
xlabel('Tempo [ms]')
legend('v_o(k)','r(k)','location','best')
ylabel('Tensão [V]')
set(gca,'fontname','times') 
set(gca,'fontsize',10)
subplot(2,1,2)
plot(rk(:,1)*1000,uk,'-k','linewidth',1); grid on;
xlim([0.053*1000,0.115*1000]);
ylim([-0.7,0.7]);
xlabel('Tempo [ms]')
legend('u(k)','location','best')
ylabel('Ação de controle - u(k)')
set(gca,'fontname','times') 
set(gca,'fontsize',10)
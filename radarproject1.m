%%%% Radar Project 
%% Ratna Sai Kiran


%% implementing P4 phased signal and calculated the autocorrelation of the P4 signal

%% generating phases
% Number of phases
M=25;

%array consisting of all phasers
M_phase=zeros(1,M);
%M_phase=mod(M_phase,2*pi)
%
for i =1:M
    M_phase(i)=((pi/M)*(i-1)*(i-1))-(pi*(i-1));
end
M_phase=mod(M_phase,2*pi)

M_phase'

%% radar paramters
MM_phase=[0 , -0.96   ,    -1.84   ,    -2.64   ,    -3.36 ,         -4  ,     -4.56 ,      -5.04  ,     -5.44   ,    -5.76  ,        -6 ,      -6.16  ,     -6.24    ,   -6.24   ,    -6.16     ,     -6      , -5.76    ,   -5.44   ,    -5.04    ,   -4.56   ,       -4 ,      -3.36 ,      -2.64    ,   -1.84    ,   -0.96]

%% speed of light  meters/sec
c_mps=3e8;

%% pulse repetition interval in Sec
PRI_sec=5;

%% carrier frequency in Hz
Fc_Hz=2e3;

%% Sampling frequency
Fs_Hz=10e4;

%% delta seconds
delta_sec=1/Fs_Hz;

%% time axis
taxis_s=0:delta_sec:PRI_sec-delta_sec;

% Number of samples
N_samp_nu=round(PRI_sec/delta_sec);


%% Inpulse duration

InPulse_sec=5;

%% number of samples when pulse is on
N_samp_inpulse=round(InPulse_sec/delta_sec);

%% Number of samples for each bit
N_samp_bits=round(N_samp_inpulse/25);
%% transmited signal 
Stx_volts=zeros(1,N_samp_nu);



%% assigning each bits its value
for i=1:M
    Stx_volts((i-1)*N_samp_bits+1:(i)*N_samp_bits)=ones(1,N_samp_bits).*exp(-j*M_phase(i));
end

Stx_volts=Stx_volts.*exp(j*2*pi*taxis_s);

% %% recieved signal
% Srx_volts=zeros(1,N_samp_nu);
% for i=1:25
%     Srx_volts(((i-1)*N_samp_bits)+1+N_samp_inpulse:((i)*N_samp_bits)+N_samp_inpulse)=ones(1,N_samp_bits).*exp(j*M_phase(i));
% end
% 
% 
% Srx_volts=Srx_volts.*exp(j*2*pi*taxis_s);

sout_volts=xcorr(Stx_volts);

% time axis after autocorrelation
t_acorr=-PRI_sec+delta_sec:delta_sec:PRI_sec-delta_sec;

pk=max(abs(sout_volts));

sout_volts_dB=20*log10(abs(sout_volts)/max(abs(sout_volts)));
t_acorr=t_acorr/0.2;
figure(1)
plot(t_acorr,(sout_volts_dB))
xlim([0,25])
ylim([-50,0])
xlabel('delay/tc')
ylabel('abs(autocorr) dB')
title('autocorrelation (dB)of P4')

figure(2)
plot(t_acorr,sout_volts_dB)
xlim([0,25])
ylim([-50,0])
xlabel('time')
ylabel('abs(autocorrelation)')
title('autocorrelation of P4')

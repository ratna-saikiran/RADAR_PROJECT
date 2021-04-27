%% Ratna Sai Kiran
%% Radar Projecyt
%% M*M MCPC Signal implementation

%%
%% M value 

M=5;

%% P4 SIGNAL PHASSE VALUES
%array consisting of all phasers
P4_phase=zeros(1,25);
%M_phase=mod(M_phase,2*pi)
%
for i =1:25
    P4_phase(i)=((pi/25)*(i-1)*(i-1))-(pi*(i-1));
end

P4_phase=mod(P4_phase,2*pi);


%% bit duration in P4 Signal
tc_sec=0.2;
%% bit duration 
tb_sec=tc_sec*M;

%%speed of light
c_mps=3e8;

%%sequences of phase values in degrees
M_phase=[0,-144,-216,-216,-144;-144,-216,-216,-144,0;-216,-216,-144,0,-144;-216,-144,0,-144,-216;-144,0,-144,-216,-216];

M_phase=M_phase.*(pi/180)
%% taken sequence was 3,5,2,4,1

%% Fs
Fc_Hz=1/tb_sec;

%% PRI 
PRI_sec=5;

%% Sampling frequency
Fs_Hz=1000;

%% delta secondas
delta_sec=1/Fs_Hz;

%% number of samples in time domain
N_samp_nu=round(PRI_sec/delta_sec);

%% number of samples where bit period
N_samp_bit=round(tb_sec/delta_sec);

%% time axis 
taxis_sec=0:delta_sec:PRI_sec-delta_sec;

%% transmitted signal
Stx_volts=zeros(1,N_samp_nu);

%% phasers order
N_phase=[3,5,2,1,4];

%% intital phase
S_phi=0;

for n=1:5
    Stx_volts(1:M*N_samp_bit)=Stx_volts(1:M*N_samp_bit)+exp(j*2*pi*Fc_Hz*(3-n)*taxis_sec(1:M*N_samp_bit));
    for m=1:M
        Stx_volts((m-1)*N_samp_bit+1:m*N_samp_bit)=Stx_volts((m-1)*N_samp_bit+1:m*N_samp_bit).*exp(j*M_phase(N_phase(n),m));
    end
end



%% time after performing auto correlation

time_acf_sec=-PRI_sec+delta_sec:delta_sec:PRI_sec-delta_sec;


Sout_volts=xcorr(Stx_volts,Stx_volts);
Sout_volts_db=10*log(abs(Sout_volts)/max(abs(Sout_volts)));

Stx_fft=fftshift(fft(Stx_volts));
%% frequency axis
F_axis=(-0.5/delta_sec)+(1/PRI_sec):1/PRI_sec:(0.5/delta_sec);


% %% fft
 figure(1)
 plot(taxis_sec,abs(Stx_volts))

%% considering only positive value
figure(2)
plot(time_acf_sec,Sout_volts_db)
xlim([0,5])
ylim([-50,0])
xlabel('delay/tb')
ylabel('abs(autocorr) dB')
title('autocorrelatoin of MCPC signal')

%% plotting the doppler ambiguity graph
%[afmg,delay,doppler]=ambgfun(Stx_volts,Fs_Hz,1/PRI_sec);
%contour(delay,doppler,afmg)
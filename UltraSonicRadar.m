
%Ultra Sonic Radar Simulation
%Rishad Ali Yasin
%YSNRIS001


%% Clear variables and command window

clear all;
close all;
clc;

%% Define constants and parameters
c = 343;                           % [m/s]  ->speed of wave

%Radar parameters
fc = 41e3;                         % [Hz]
f0 = 40e3;
T =100e-3;                           % Pulse length in [s]  
fs = 192e3;                       % [Hz]   ->Sampling rate. So 44100 samples are obtained in one second
ts = 1/fs;                         % Sampling period
B = 1e3;                           % Bandwidth (Hz) 
PRI = 200e-3;                           % seconds
NumPulses = 16;                    % Number of pulses 
R_max = c*PRI/2;                   % Maximum range of target to simulate
 
PRF = 1/PRI;                       % Hz
t_max = PRI*NumPulses;             % Maximum time to simulate 
K = B/T;                           % Chirp rate
t  = (0:ts:(t_max-ts));            % Time vector (in seconds)
NumSamplesInPRI = round(PRI*fs);
 
% Target parameters 
R_target = 1.4;                       % Target range 
Vel_target = 0.005;                  % Target radial velocity in m/s
td = 2* R_target/c;                  % Time taken for signal to travel to target and back to the radar
 
 
% Simulation parameters
InitialDelay_s = 1;                   % seconds
clims = [-50 0 ];

%% Calculate and display radar parameters
lamda = c/fc;
PulseCompressionGain = T*B;
BlindRange_m = c*T/2/PulseCompressionGain;
UnambiguousRange_m = c*PRI/2;
RangeResolution_m = c/(2*B);
UnambiguousVelocity_cms = PRF/2*lamda/2*100;
 
clc
disp(' ');
disp(['Blind Range = ' num2str(roundn(BlindRange_m,-3)) ' m']);
disp(['Unambiguous Range = ' num2str(roundn(UnambiguousRange_m,0)) ' m']);
disp(['Range Resolution = ' num2str(roundn(RangeResolution_m,-3)) ' m']);
disp(['Unambiguous Velocity = ' num2str(roundn(UnambiguousVelocity_cms,-2)) ' cm/s']);
disp(['Pulse Compression gain = ' num2str(roundn(PulseCompressionGain,0)) '']);
disp(' ');

%% Generate the transmit pulse
% Generate Transmit signal 
Tx_Signal = cos(2*pi*(f0*(t - T/2) + 0.5*K*(t - T/2).^2) ).* rect( (t - T/2)/T )   ;  % transmit signal with chirp pulse
    
% Generate the transmit pulse 
NumSamplesTxPulse = ceil(T/ts);             % number of samples of the transmit pulse 
Tx_p = Tx_Signal(1: NumSamplesTxPulse);     % transmit pulse only
FreqAxis = (-NumSamplesTxPulse/2:1:NumSamplesTxPulse/2-1)*fs/NumSamplesTxPulse;
t_Tx_p = t(1: NumSamplesTxPulse);
figure(1);
plot(t_Tx_p,Tx_p);
figure(2);
plot(FreqAxis,20*log10(fftshift(abs(fft(Tx_p)))));


%% Generate the range line for multiple pulses, assuming a stationary target
 Rx_Signal = zeros(1, size(t,2));
 
 for Count_PulseNum = 1: NumPulses
 
    tdn = 2*(R_target - Vel_target*PRI*(Count_PulseNum-1))/c + PRI*(Count_PulseNum-1); 
     
    Rx_Signal = Rx_Signal +  cos(2*pi*(f0*(t - T/2 - tdn) + 0.5*K*(t - T/2 - tdn).^2) ).* rect( (t - T/2 - tdn)/T );  % received signal
    
    %figure;
    %plot(Rx_Signal);
 
 end
N2 = size(Rx_Signal,2); %Number of samples in the received signal
t2 = (0:1:(N2-1))*ts; %Time vecctor for received signal
figure();
plot(t2,Rx_Signal)
FreqAxis = (-size(Rx_Signal,2)/2:1:size(Rx_Signal,2)/2-1)*fs/size(Rx_Signal,2);
figure();
plot(FreqAxis,20*log10(fftshift(abs(fft(Rx_Signal)))));


%% downmix both the transmit pulse and the received signal
t_Tx_p = t(1: NumSamplesTxPulse);
I_channel = Tx_p.*cos(2*pi*fc*(t_Tx_p));
Q_channel = Tx_p.*-sin(2*pi*fc*(t_Tx_p));
y = I_channel + 1i*Q_channel;
load('UltraSonic_coefficients.mat');
Signal_I_AferLPF = filter(b,1,I_channel);
Signal_Q_AferLPF = filter(b,1,Q_channel);
y1 = complex(Signal_I_AferLPF,Signal_Q_AferLPF);%Base band of transmit pulse
FreqAxis = (-NumSamplesTxPulse/2:1:NumSamplesTxPulse/2-1)*fs/NumSamplesTxPulse;
figure;
plot(FreqAxis,fftshift(20*log10(abs(fft(y1)))))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_channel1 =  Rx_Signal.*cos(2*pi*fc*t2);
Q_channel1 = Rx_Signal.*-sin(2*pi*fc*t2); 
z = I_channel1 + 1i*Q_channel1;
Signal_I_AferLPF1 = filter(b,1,I_channel1);
Signal_Q_AferLPF1 = filter(b,1,Q_channel1);
RxSignal_filtered = complex(Signal_I_AferLPF1,Signal_Q_AferLPF1);
figure;
FreqAxis1 = (-N2/2:1:N2/2-1)*fs/N2;
plot(FreqAxis1,fftshift(20*log10(abs(fft(RxSignal_filtered)))))

%% Re-shape receive signal into a matrix
RxSignalMatrix1 = reshape(RxSignal_filtered,  length(RxSignal_filtered)/NumPulses, NumPulses);
RxSignalMatrix = RxSignalMatrix1.'; % Each row is a range line


%% Pulse compression or matched filtering

%Do matched filtering in the time-domain
%y_time1 = conv(h,RxSignal_filtered);
%figure; plot(abs(y_time1));

% Frequency domain: matched filtering

h = conj(fliplr([ y1  zeros(1,(size(RxSignalMatrix,2)- NumSamplesTxPulse))])); %Matched filter
h1 = h;
H = fft(h1); %Take the FFT of the matched filer
H = repmat(H,NumPulses,1);%Replicte the matched filter into a column vector
%RxSignal = transpose(y);
%figure;
%imagesc(20*log10(abs(RxSignal_r)));
%figure;
Z = fft(RxSignalMatrix,[],2); %Take FFT along every row of the RxSignal_r matrix
%plot(20*log10(abs(RxSignal_r(1,:))));
Y = H.*Z; %Multiply in frequency domain
y_time = ifft(Y,[],2); %Take inverse FFT along evey row
M = size(RxSignalMatrix,2); %Number of columns in the RxSignal_r matrix
Range = (0:1:(M-1))*ts*c/2; %Range vector
Num_Pulses = 1:1:NumPulses; %Number of Pulses axis
figure;
plot(Range,20*log10(abs(y_time(1,:))));
 
  % Plot range line after pulse compression

figure;
imagesc(Range, Num_Pulses,20*log10(abs(y_time)));

  
  %% Generate the Range-Doppler map
  
  % Obtain the window matrix (A) and multiply element-by-element with the RangeLines (B)
  %         C = A.*B
  % Then apply the FFT in the slow-time dimension 
  %         D = fft(C, [], 1) % Each row of C is a range-line 
  
  
  % Apply the FFT in the slow-time dimension 

  
  % Plot Range-Doppler matrix
w = hamming(NumPulses);
B = repmat(w,1,M);
D = y_time.*B;
Range_Doppler = fftshift(fft(D,[],1), 1);
Range_Doppler_dB = NormLimitDynRange( Range_Doppler,40);
Doppler_freq = (-NumPulses/2:1:(NumPulses/2-1))*PRF/NumPulses; %Doppler frequency vector
velocity = Doppler_freq*lamda/2; %Velocity axis
figure;
imagesc(Range,velocity,(Range_Doppler_dB));
colormap('jet');  



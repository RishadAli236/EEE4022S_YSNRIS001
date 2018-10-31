%Ultrasonic Radar Range Doppler Processing
%Rishad Ali Yasin
%YSNRIS001

%% Clear variables and command window

clear all;
close all;
clc;

%% Define constants and parameters

c = 343;                           % Speed of sound [m/s]
f0 = 40e3;                         % Centre frequency [Hz]
T = 100e-3;                         % Pulse width [s]  
fs = 192e3;                        % Sampling rate [Hz]
B = 1e3;                             % Bandwidth [Hz] 
PRI = 200e-3;                       % Pulse Repetition Interval [s]
NumPulses = 16;                    % Number of pulses 
R_max = c*PRI/2;                   % Unambiguous range [m]
PRF = 1/PRI;                       % Pulse Repetition Frequency [Hz]
wave_length = c/f0;                %Wave length

%% Generate Transmit Pulse and Transmit

Tx_pulse = pulse_gen(T,fs,f0,B);
Tx_pulse = Tx_pulse - mean(Tx_pulse); 

%% Record Echo

recObj = audiorecorder(fs,16,1);
obj = serial('COM4','BaudRate',12000000,'DataBits',8);
fopen(obj);
fprintf(obj,'1');
recordblocking(recObj, PRI*(NumPulses+4));  % Records audio for a fixed number of seconds
fclose(obj);
delete(obj);
Rx_sig = transpose(getaudiodata(recObj,'double'));

Rx_sig = Rx_sig(1:PRI*fs*16);



ts = 1/fs;
t_Rx = (0:ts:(size(Rx_sig,2)*ts)-ts);
FreqAxis = (-size(Rx_sig,2)/2:1:size(Rx_sig,2)/2-1)*fs/size(Rx_sig,2);
figure();
plot(t_Rx,Rx_sig);
figure();
plot(FreqAxis,20*log10(fftshift(abs(fft(Rx_sig)))));

%% Band pass filtering
load('Bandpass.mat');
Rx_sig = filter(Band_pass,1,Rx_sig); %Apply band pass filter to received signal
figure();
plot(FreqAxis,20*log10(fftshift(abs(fft(Rx_sig)))));
% 
% 
%% Down Mixing

t_Tx = (0:ts:T-ts);
I_channel = Tx_pulse.*cos(2*pi*f0*(t_Tx)); %I-channel
Q_channel = Tx_pulse.*-sin(2*pi*f0*(t_Tx)); %Q-channel
y = I_channel + 1i*Q_channel;
load('UltraSonic_coefficients.mat'); % Low pass filter
Signal_I_AferLPF = filter(b,1,I_channel);
Signal_Q_AferLPF = filter(b,1,Q_channel);
Tx_pulse_baseband = complex(Signal_I_AferLPF,Signal_Q_AferLPF);%Base band of transmit pulse
FreqAxis = (-size(Tx_pulse,2)/2:1:size(Tx_pulse,2)/2-1)*fs/size(Tx_pulse,2);
figure();
plot(FreqAxis,fftshift(20*log10(abs(fft(Tx_pulse_baseband)))))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_channel1 =  Rx_sig.*cos(2*pi*f0*t_Rx);
Q_channel1 = Rx_sig.*-sin(2*pi*f0*t_Rx); 
z = I_channel1 + 1i*Q_channel1;
Signal_I_AferLPF1 = filter(b,1,I_channel1);
Signal_Q_AferLPF1 = filter(b,1,Q_channel1);
Rx_sig_baseband = complex(Signal_I_AferLPF1,Signal_Q_AferLPF1);
figure();
FreqAxis = (-size(Rx_sig,2)/2:1:size(Rx_sig,2)/2-1)*fs/size(Rx_sig,2);
plot(FreqAxis,fftshift(20*log10(abs(fft(Rx_sig_baseband)))));
% 
% 
% %% Re-shape receive signal into a matrix
 RxSignalMatrix = transpose(reshape(Rx_sig_baseband,  length(Rx_sig_baseband)/NumPulses, NumPulses));
% 
% 
 %% Exercise 4/6: Pulse compression or matched filtering

 h = conj(fliplr([ Tx_pulse_baseband  zeros(1,(size(RxSignalMatrix,2)- size(Tx_pulse_baseband,2)))])); %Matched filter Padded Matched filter with zeros so 
                                                                                                        %that it is the same size as RxSignal_r
% Frequency domain: matched filtering
h1 = h;
H = fft(h1); %Take the FFT of the matched filer
H = repmat(H,NumPulses,1);%Replicte the matched filter into a column vector
Z = fft(RxSignalMatrix,[],2); %Take FFT along every row of the RxSignal_r matrix
Y = H.*Z; %Multiply in frequency domain
y_time = ifft(Y,[],2); %Take inverse FFT along evey row
 M = size(RxSignalMatrix,2); %Number of columns in the RxSignal_r matrix
 Range = (0:1:(M-1))*ts*c/2; %Range vector
 figure(8);
 plot(Range,(abs(y_time(1,:))));
Num_Pulses = 1:1:NumPulses; %Number of Pulses axis
%  
%   
% %   % Plot range line after pulse compression
% % 


% Code to do compensation in phase. Make leakage appear at Velocity = 0 m/s

AvgRangeLine = abs(y_time(1,:));
%[MaxVal, MaxIndx] = max(AvgRangeLine);
MaxIndx = find(AvgRangeLine>10);
MaxIndx = MaxIndx(1);
%figure; plot(AvgRangeLine); grid on;

Leakage = y_time(:, MaxIndx);
CompensationVec = conj(Leakage);
CompensationMatrix = repmat(CompensationVec, 1, size(y_time,2)); 

y_time_AfterComp = y_time.*CompensationMatrix;

y_time = circshift(y_time_AfterComp,-MaxIndx, 2);
figure();
%imagesc(Range(1000:(7000)), Num_Pulses,((abs(y_time(:, 1000:(7000))))));
imagesc(Range, Num_Pulses,((abs(y_time))));
colorbar; 
colormap('jet');

% End of code for compensation

%% Genereate Range Doppler Map
 w = hamming(NumPulses); %Apply hamming windowing function
 B = repmat(w,1,M);
 D = y_time.*B;
 Range_Doppler = fftshift(fft(D,[],1), 1); %Apply FFT along every column 
 Doppler_freq = (-NumPulses/2:1:(NumPulses/2-1))*PRF/NumPulses; %Doppler frequency vector
 velocity = Doppler_freq*wave_length/2; %Velocity axis
 figure;
 imagesc(Range(1000:(7000)),velocity,(abs(Range_Doppler(:, 1000:(7000)))));axis ij;
 xlabel('Range (m)', 'fontsize', 12);
 ylabel('Velocity (m/s)', 'fontsize', 12);
 colorbar;
 colormap('jet');


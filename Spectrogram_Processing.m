%Ultrasonic Radar Spectrogram Proceessing
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
figure(2);
plot(t_Rx,Rx_sig);
figure(3);
plot(FreqAxis,20*log10(fftshift(abs(fft(Rx_sig)))));

%% Generate Sectrogram
w0 = 4.04e4/(192e3/2);
BW = 240/(192e3/2);
[b,a] = iirnotch(w0,BW);  %Generate a notch filter
Rx_sig = filter(b,a,Rx_sig); %Apply the notch filter
[S, f, t, P] = spectrogram(Rx_sig, 2^13, 2^13-50, [], 192e3, 'yaxis');

%% Limit Range of Doppler Frequencies and Dynamic Range 
f1 = f > 38e3 & f < 42e3;
P1 = P(f1, :);
f1 = f(f1);
 P1 = P1./(max(max(P1))); 
 clims = [-80  0];
figure; imagesc(t, f1, 20*log10((P1)), clims);
colormap('jet');
colormap;
axis ij; 





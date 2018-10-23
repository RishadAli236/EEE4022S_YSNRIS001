function [ y ] = pulse_gen( T,fs,f0,B )
ts = 1/fs;
K = B/T;
bits = 12;
t = (0:ts:T-ts);
%y = sin (2*pi*f0*t);
y = cos(2*pi*(f0*(t - T/2) + 0.5*K*(t - T/2).^2) ).* rect( (t - T/2)/T )   ;
%y = chirp(t, f0, T-ts, f0+B);
% y = y + 1;
% y = y*(((2^bits)-1)/2);
% y = round(y); 
figure();
plot(t,y); grid
FreqAxis = (-size(y,2)/2:1:size(y,2)/2-1)*fs/size(y,2);
figure();
plot(FreqAxis,20*log10(fftshift(abs(fft(y)))));
end


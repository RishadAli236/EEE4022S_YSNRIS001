% Code to strem data from ADC
% Rishad Ali Yasin
% YSNRIS001


%%
obj = serial('COM4','BaudRate',12000000,'DataBits',8);
obj.InputBufferSize = 5000000; 
obj.Timeout = 10;
fopen(obj);
fprintf(obj,'1');
x = fread(obj,4486,'uint8');
x = transpose(x);
disp(obj)
fclose(obj);
delete(obj);
Rx_sig = zeros(1,2243);
i = 1;
for n = 1:2243
   Rx_sig(n) = bitor(bitshift(x(i),8),x(i+1)); 
   i = i + 2;
end
fs = 89743;
ts = 1/fs;
t_Rx = (0:ts:(size(Rx_sig,2)*ts)-ts);
FreqAxis = (-size(Rx_sig,2)/2:1:size(Rx_sig,2)/2-1)*fs/size(Rx_sig,2);
figure();
plot(t_Rx,Rx_sig); grid
figure();
plot(FreqAxis,20*log10(fftshift(abs(fft(Rx_sig)))));


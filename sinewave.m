%% Generating an n-bit sine wave
% Modifiable parameters: step, bits, offset
close; clear; clc;

points = 32;                            % number of points between sin(0) to sin(2*pi)
bits   = 12;                             % 12-bit sine wave for 12-bit DAC


t = (0:((2*pi/(points-1))):(2*pi))*1;        % creating a vector from 0 to 2*pi
y = sin(t);                              % getting the sine values
y = y + 1;                               % getting rid of negative values (shifting up by 1)
y = y*(((2^bits)-1)/2);
y = round(y);                            % rounding the values
plot(t, y); grid                         % plotting for visual confirmation

fprintf('%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, \n', y);
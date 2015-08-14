clear
Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = 1000;                     % Length of signal
t = (0:L-1)*T;                % Time vector
% Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
x = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t); 
y = x %+ 2*randn(size(t));     % Sinusoids plus noise
%plot(Fs*t(1:50),y(1:50))
%plot(t,y)
Y=fft(y);
plot((0:501),abs(Y(1:502)));

s=ifft(Y);
plot(t(1:50),s(1:50))
hold on
plot(t(1:50),y(1:50)+0.1,'red')



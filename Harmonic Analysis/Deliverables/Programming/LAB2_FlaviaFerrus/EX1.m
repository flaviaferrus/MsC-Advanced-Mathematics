%% LAB 2: The Z-transform and filter design
% Flàvia Ferrús Marimón

%% EXAMPLE 1

% Let's consider an initial FIR filter 
% by windowing using the function fir1

% Low pass band with 48º and 
% 0.35 pi < omega < 0.65pi rad

b = fir1(48, [0.35 0.65]);
freqz(b,1,512)

%c = fir2(48, [0.35 0.65]);
%freqz(c,1,512)


% Consider a custom low pass filter 
% on a generated signal:
load chirp.mat

t = (0:length(y) - 1)/Fs;
blo = fir1(34, 0.48, chebwin(35, 30));

outlo = filter(blo, 1, y);
subplot(2,1,1)
plot(t, y)
title('Original Signal')
ys = ylim;

subplot(2,1,2)
plot(t, outlo)
title('Lowpass Filtered Signal')
xlabel('Time (s)')
ylim(ys)

blo = fir1(34, 0.48, chebwin(35, 30));
outlo_b = filter(b, 1, y);
subplot(2,1,1)
plot(t, y)
title('Original Signal')

ys = ylim;
subplot(2,1,2)
plot(t, outlo_b)
title('Lowpass Filtered Signal')
xlabel('Time (s)')
ylim(ys)

freqz(blo,1,512)


%% EXERCISE 1 
% Not used: 

%% Using the fir1 function:

% Define the values of theta and the filter order
theta = pi/6;
N = 2;
theta_d = 180/6;
% Compute and plot the frequency response using freqz
%[H, omega] = freqz(b, 1, 1000, 'whole');
%magnitude_H = abs(H);

% Plot the magnitude of the frequency response
%plot(omega, magnitude_H);
%title('Magnitude of the Filter Frequency Response');
%xlabel('\omega');
%ylabel('|H(e^{i\omega})|');
%axis([-pi pi 0 1.5]);
%grid on;


blo = fir1(theta, [0.35 0.65]);
[H, omega] = freqz(blo,1,512);
plot(omega, abs(H));
title('Magnitude of the Filter Frequency Response');
xlabel('\omega');
ylabel('|H(e^{i\omega})|');
axis([-pi pi 0 1.5]);
grid on;



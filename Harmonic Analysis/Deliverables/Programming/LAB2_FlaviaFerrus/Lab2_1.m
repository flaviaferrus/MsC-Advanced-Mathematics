% Define the values of theta and omega
theta = pi/6;
omega = linspace(-pi, pi, 1000);

% Compute the frequency response
H = 1 - 2*cos(theta)*exp(-1i*omega) + exp(-2i*omega);
magnitude_H = abs(H);

% Plot the magnitude of the frequency response
plot(omega, magnitude_H);
title('Magnitude of the Filter Frequency Response');
xlabel('\omega');
ylabel('|H(e^{i\omega})|');
axis([-pi pi 0 1.5]);
grid on;

figure
freqz(filter_h(theta))


%% LOADING AUDIO SIGNALS

[x,FS] = audioread('easy.wav');
% Hear the original signal
player = audioplayer(0.8*x, FS);
play (player);

figure 
plot(1: length(0.8*x), 0.8*x);
title('Original Signal');
xlabel('Hz');
ylabel('Signal');


% Add sinusoidal noise to the signal
t = (1:length(x))/FS;
y = 0.8*x + 0.1*sin(35000* t');
% Listen to the signal again
player = audioplayer (y, FS);
play (player);

figure 
plot(1: length(y), y);
title('Signal with noise');
xlabel('Hz');
ylabel('Signal');


%% Filter using lowpass filter

% Filter the noisy signal

y_filtered = lowpass(y,5,FS);
y_filtered2 = lowpass(y, 350, FS, ImpulseResponse= "fir");
y_filtered3 = lowpass(y, 15000, FS, ImpulseResponse= "fir");

% Play the filtered signal
player = audioplayer(y_filtered2, FS);
play(player);


figure 
plot(1: length(y_filtered2), y_filtered2);
title('Filtered Signal');
xlabel('Hz');
ylabel('Signal');

figure 
plot(1: length(y_filtered3), y_filtered3);
title('Filtered Signal');
xlabel('Hz');
ylabel('Signal');


%% Filter using H

% Design the filter with two zeros
theta_ = (35000) / FS;
%theta_ = pi/6; %2*pi*35000
fprintf("Initial frequency FS= %.2f\n", FS);
fprintf("Value of theta= %.2f\n", theta_);

h_filter = [1, -2*cos(theta_), 1];

% Filter the noisy signal
y_filtered_H = filter(h_filter, 1, y);

% Play the filtered signal
player = audioplayer(y_filtered_H, FS);
play(player);

figure 
plot(1: length(y_filtered_H), y_filtered_H);
title('Filtered Signal using H');
xlabel('Hz');
ylabel('Signal');

y_filtered_H2 = conv(y, filter_h(theta_));
figure 
plot(1: length(y_filtered_H2), y_filtered_H2);
title('Filtered Signal with H');
xlabel('Hz');
ylabel('Signal');



% Function to compute filter
function h_ = filter_h(theta_)
    h_ = [1, -2*cos(theta_), 1];
end


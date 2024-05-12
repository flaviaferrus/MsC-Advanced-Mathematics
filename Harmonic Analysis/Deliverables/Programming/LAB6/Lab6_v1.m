[x,FS] = audioread('easy.wav');
player = audioplayer (0.8*x, FS);
play (player);

y = x + normrnd(0, 0.2, length(x),1);
player = audioplayer(y, FS);
play (player);

% Ratio to keep
r = 0.9;
% Perform the discrete wavelet transform
wa = 'db6';
J = 5;
[c, S] = wavedec(y, J, wa);

% Calculate the number of coefficients to keep based on the ratio
numCoeffsToKeep = round(r * numel(c));

% Sort the coefficients in descending order of magnitude
[~, sortedIndices] = sort(abs(c), 'descend');

% Keep the most significant coefficients
cc = zeros(size(c));
cc(sortedIndices(1:numCoeffsToKeep)) = c(sortedIndices(1:numCoeffsToKeep));

% Reconstruct the audio signal using the inverse wavelet transform
y_reconstructed = waverec(cc, S, wa);

% Play the attenuated audio signal
player = audioplayer(y_reconstructed, FS);
play(player);

figure 
plot(1: length(y_reconstructed), y_reconstructed);
title('Filtered Signal with H');
xlabel('Hz');
ylabel('Signal');

figure 
plot(1: length(y), y);
title('Noisy Signal');
xlabel('Hz');
ylabel('Signal');

figure 
plot(1: length(x), x);
title('Original Signal');
xlabel('Hz');
ylabel('Signal');


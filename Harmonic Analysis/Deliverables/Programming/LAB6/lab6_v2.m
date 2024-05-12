% Add the LTFAT toolbox path
addpath('ltfat');

% Read the cameraman image
f = imread('cameraman.tif');

% Ratio to keep
r = 0.05;

%% Parameters for the Wavelet system
% Analysis filters
wa = 'db6';
% Synthesis filters
ws = 'db6';
% No. of levels
J = 5;

%% Show the original image
figure(1);
imagesc(f);
colormap(gray);
axis('image');
title('Original Image');

%% Compressed images
figure(2);

% Perform the wavelet transform using 'wavedec2'
[c_fwt, S] = wavedec2(f, J, wa);

% Calculate the number of coefficients to keep based on the ratio
numCoeffsToKeep = round(r * numel(c_fwt));

% Sort the coefficients in descending order of magnitude
[~, sortedIndices] = sort(abs(c_fwt), 'descend');

% Keep the largest coefficients
cc_fwt = zeros(size(c_fwt));
cc_fwt(sortedIndices(1:numCoeffsToKeep)) = c_fwt(sortedIndices(1:numCoeffsToKeep));

% Reconstruct the image using the inverse wavelet transform 'waverec2'
r_fwt = waverec2(cc_fwt, S, ws);

subplot(1, 2, 1);
imagesc(r_fwt);
colormap(gray);
axis('image');
title('Wavelet Compression');

% Perform the DCT on the original image
c_dct = dct2(f);

% Keep the largest coefficients based on the ratio
numCoeffsToKeep_dct = round(r * numel(c_dct));
[~, sortedIndices_dct] = sort(abs(c_dct(:)), 'descend');
cc_dct = zeros(size(c_dct));
cc_dct(sortedIndices_dct(1:numCoeffsToKeep_dct)) = c_dct(sortedIndices_dct(1:numCoeffsToKeep_dct));

% Reconstruct the image using the inverse DCT 'idct2'
r_dct = idct2(cc_dct);

subplot(1, 2, 2);
imagesc(r_dct);
colormap(gray);
axis('image');
title('DCT Compression');

% Remove the LTFAT toolbox path
rmpath('ltfat');

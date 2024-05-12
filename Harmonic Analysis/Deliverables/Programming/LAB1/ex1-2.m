% Lab 1: Basic Digital Signal Processing
% Author: Àlex Pujol Vidal





% ===> EXERCISE 1: Filtering over COVID-19 data <===

% File path and name
file_path = 'newcases.txt';


%---part1
%---Filter: Moving average

% Create figure
figure('Position', [100, 200, 600, 800]);

% Load the content of the file_in_loadpath
array_data = load(file_path);

% Plot original data
subplot(3,1,1);
x_values = 1:length(array_data);
plot(x_values, array_data);

% Add labels and title
xlabel('Days');
ylabel('New cases');
title('Evolution of new COVID-19 cases');


% Define moving average filter using for loop
function result = h_for(array_data, n)
  len = length(array_data);
  result = array_data(1:len-(n-1));
  for i = 1:(n-1)
    result = result + array_data(1+i:len-(n-1-i));
  endfor
  result /= 5;
endfunction

% Filter data with for
array_filtered = h_for(array_data, 5);

% Plot filtered data
subplot(3,1,2);
x_values = 1:length(array_filtered);
plot(x_values, array_filtered);

% Add labels and title
xlabel('Days');
ylabel('New cases');
% title('Moving average evolution of new COVID cases');


% Define a function to use as a filter for convolution
function result = h(n)
  result = ones(1,n)/n;
endfunction

% Make the convolution
array_filtered = conv(array_data, h(5));

% Plot filtered data
subplot(3,1,3);
x_values = 1:length(array_filtered);
plot(x_values, array_filtered);

% Add labels and title
xlabel('Days');
ylabel('New cases');
% title('Moving average evolution of new COVID cases');


%---part2
%---Compute and draw the Fourier transform of h

% Create figure
figure('Position', [100, 100, 400, 400]);
freqz(h(5))


%---part3
%---Alternative filters proposal

% Create figure 
figure('Position', [100, 200, 600, 1200]);

% Plot original data
subplot(5,1,1);
x_values = 1:length(array_data);
plot(x_values, array_data);

% Add labels and title
xlabel('Days');
ylabel('New cases');
title('Evolution of new COVID cases. Moving average.');

% Plot filtered data order 3 moving average
array_filtered = conv(array_data, h(3));
subplot(5,1,2);
x_values = 1:length(array_filtered);
plot(x_values, array_filtered);
xlabel('Days');
ylabel('New cases');

% Plot filtered data order 5 moving average
array_filtered = conv(array_data, h(5));
subplot(5,1,3);
x_values = 1:length(array_filtered);
plot(x_values, array_filtered);
xlabel('Days');
ylabel('New cases');

% Plot filtered data order 7 moving average
array_filtered = conv(array_data, h(7));
subplot(5,1,4);
x_values = 1:length(array_filtered);
plot(x_values, array_filtered);
xlabel('Days');
ylabel('New cases');

% Plot filtered data order 11 moving average
array_filtered = conv(array_data, h(11));
subplot(5,1,5);
x_values = 1:length(array_filtered);
plot(x_values, array_filtered);
xlabel('Days');
ylabel('New cases');



%---Other filters
%---Gaussian filter

% Create figure 
figure('Position', [100, 200, 600, 1200]);

function smoothed_data = gaussian_filter(data, sigma)
  % data: input time series data
  % sigma: standard deviation of the Gaussian distribution

  % Compute the size of the filter window
  window_size = ceil(6 * sigma);
  
  % Create the filter window centered at zero
  filter_window = -floor(window_size/2):floor(window_size/2);

  % Compute the Gaussian filter
  filter = exp(-(filter_window.^2) / (2 * sigma^2));
  filter = filter / sum(filter);

  % Apply the filter to the data using convolution
  smoothed_data = conv(data, filter, 'same');
end

% Plot original data
subplot(5,1,1);
x_values = 1:length(array_data);
plot(x_values, array_data);

% Add labels and title
xlabel('Days');
ylabel('New cases');
title('Evolution of new COVID cases. Gaussian filter.');

% Plot filtered data of Gaussian sigma 0.5
array_filtered = gaussian_filter(array_data, 0.5);
subplot(5,1,2);
x_values = 1:length(array_filtered);
plot(x_values, array_filtered);
xlabel('Days');
ylabel('New cases');

% Plot filtered data of Gaussian sigma 1
array_filtered = gaussian_filter(array_data, 1);
subplot(5,1,3);
x_values = 1:length(array_filtered);
plot(x_values, array_filtered);
xlabel('Days');
ylabel('New cases');

% Plot filtered data of Gaussian sigma 1.5
array_filtered = gaussian_filter(array_data, 1.5);
subplot(5,1,4);
x_values = 1:length(array_filtered);
plot(x_values, array_filtered);
xlabel('Days');
ylabel('New cases');

% Plot filtered data of Gaussian sigma 2
array_filtered = gaussian_filter(array_data, 2);
subplot(5,1,5);
x_values = 1:length(array_filtered);
plot(x_values, array_filtered);
xlabel('Days');
ylabel('New cases');



% ===> EXERCISE 2: Removing noise from Wikiloc's GPS signal data <===

% File path and name
file_path = 'heights.txt';

% Create figure
figure('Position', [100, 200, 600, 800]);

% Load the content of the file_in_loadpath
array_data = load(file_path);

% Plot original data
subplot(3,1,1);
x_values = 1:length(array_data);
plot(x_values, array_data);

% Add labels and title
xlabel('Distance');
ylabel('Height');
title('GPS height signal');

% Compute climbing and descending
function [climb, desc] = slope(x)
  len = length(x);
  climb = 0;
  desc = 0;
  for i = 1:(len-1)
    if x(i+1)>x(i)
      climb = climb + x(i+1) - x(i);
    else
      desc = desc + x(i) - x(i+1);
    end      
  endfor
endfunction

[climb, desc] = slope(array_data);

output_str = sprintf("Original data:\nClimbing = %.2f\nDescending = %.2f\n", climb, desc);
disp(output_str)

% Applying moving average of order 5
array_filtered = conv(array_data, h(5));

% Recomputing climbing and descending
[climb, desc] = slope(array_filtered(5:length(array_filtered)-5));

output_str = sprintf("After Moving Average:\nClimbing = %.2f\nDescending = %.2f\n", climb, desc);
disp(output_str)

% Plotting filtered data
subplot(3,1,2);
x_values = 5:length(array_filtered)-5;
plot(x_values, array_filtered(5:length(array_filtered)-5));
xlabel('Distance');
ylabel('Height');

% Applying Gaussian filter with sigma 2
array_filtered = gaussian_filter(array_data, 2);

% Recomputing climbing and descending
[climb, desc] = slope(array_filtered(5:length(array_filtered)-5));

output_str = sprintf("After Gaussian filter:\nClimbing = %.2f\nDescending = %.2f\n", climb, desc);
disp(output_str)

% Plotting filtered data
subplot(3,1,3);
x_values = 5:length(array_filtered)-5;
plot(x_values, array_filtered(5:length(array_filtered)-5));
xlabel('Distance');
ylabel('Height');

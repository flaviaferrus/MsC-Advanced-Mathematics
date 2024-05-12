%% LAB 3: FFT

%% EXERCISE 1
k = [0,1,2,3,4,5,6,7];
fprintf("%.2f\n", k);

% MATLAB indexes start at 1 not 0
% Then the even are actually the odd (we start with a_0 at 
% index position 1

fprintf("Length : %.2f\n", length(k));

fprintf("----Computing the RFFT-----\n");
fprintf("Custom function:\n");
k_RFFT = RFFT(k);
fprintf("%f%+fj\n", real(k_RFFT), imag(k_RFFT));
fprintf("Matlab function:\n")
k_fft_m = fft(k);
fprintf("%f%+fj\n", real(k_fft_m), imag(k_fft_m));


fprintf("----Computing the IFFT-----\n");
fprintf("Custom function:\n");
k_IFFT = RIFFT__(k_RFFT)/length(k_RFFT);
%fprintf("%f%+fj\n", real(k_IFFT), imag(k_IFFT));
fprintf("%.2f\n", k_IFFT);

fprintf("Custom function no recursive:\n");
k_IFFT_2 = IFFT__(k_RFFT);
%fprintf("%f%+fj\n", real(k_IFFT_2), imag(k_IFFT_2));
fprintf("%.2f\n", k_IFFT_2);

fprintf("Matlab function:\n")
k_ifft_m = ifft(k_RFFT);
%fprintf("%f%+fj\n", real(k_ifft_m), imag(k_ifft_m));
fprintf("%.2f\n", k_ifft_m);

fprintf("Matlab function 2:\n")
k_ifft_m2 = ifft(k_fft_m);
%fprintf("%f%+fj\n", real(k_ifft_m2), imag(k_ifft_m2));
fprintf("%.2f\n", k_ifft_m2);


function IFFT_ = RIFFT__(a)
    n = length(a);
    if n == 1 
        IFFT_ = a;
        return
    end 
    w_n = exp(2*pi*1i/n); % Note the positive sign here
    w = 1;
    a_even = a(1:2:end);
    a_odd = a(2:2:end);
    y_even = RIFFT__(a_even);
    y_odd = RIFFT__(a_odd);
    y = zeros(1, n);
    for k = 1:(floor(n/2) )
        y(k) = y_even(k) + w*y_odd(k);
        y(k + floor(n/2)) = y_even(k) - w*y_odd(k);
        w = w*w_n;
    end 
    IFFT_ = y; % Divide by n to incorporate the normalization factor
end

function x = IFFT__(X)
    X = conj(X);
    n = length(X);
    x = RFFT(X)/n;
    x = conj(x);
end


function RFFT_ = RFFT(a)
    n = length(a);
    if n == 1 
        RFFT_ = a;
        return
    end 
    w_n = exp(-2*pi*1i/n);
    w = 1;
    a_even = a(1:2:(n));
    a_odd = a(2:2:(n));
    y_even = RFFT(a_even);
    y_odd = RFFT(a_odd);
    y = zeros(1, n);
    for k = 1:(floor(n/2))
        y(k) = y_even(k) + w*y_odd(k);
        y(k + floor(n/2)) = y_even(k) - w*y_odd(k);
        w = w*w_n;
    end 
    RFFT_ = y;
end


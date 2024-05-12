%% LAB 3: FFT

%% EXERCISE 2

N = 8;

suma = 0;
i = 2;
p = [];
%p = zeros(1, 2*N);
while suma < N 
    if isprime(i)
        p(i) = 1;
        fprintf("%.2f\n", i);
    else 
        p(i) = 0;
    end
    suma = sum(p);
    i = i+1;
end

fprintf("Polynomial p: \n");
fprintf("%.2f\n", p);

%% Now we have a polynomial, of degree M, with the first N
% primes. We seek to consider the point-value representation
% Therefore, we may consider the polynomial of degree 2M
% by adding the M zero coefficients to p


M = length(p);
p_2 = zeros(1, M);
p_3 = cat(2, p, p_2);
fprintf("Polynomial p with zeros: \n");
fprintf("%.2f\n", p_3);

%p_fft = RFFT(p_3);
p_fft = fft(p_3);
q_fft = p_fft .* p_fft;
fprintf("Polynomial q: \n");
fprintf("%.2f\n", q_fft);

%q_ = RIFFT(q_fft);
q_ = ifft(q_fft);

fprintf("Polynomial q: \n");
fprintf("%.2f\n", q_);
%fprintf("%.2f\n", real(q_));
%fprintf("%.2f\n", imag(q_));
%disp(length(q_));
disp(q_)
q_2 = int64(q_);
q_trip_3 = q_2(mod(q_2, 2) == 1);
%fprintf("Odd numbers from q_:\n");
%fprintf("%.2f\n",mod(q_2, 2) == 1);

q_3 = q_2 == 3;
fprintf("3's q_:\n");
%fprintf("%.2f\n", q_2(q_3));
%fprintf("%.2f\n", q_3);

q_trip_3_ind = q_trip_3 >= 3;
fprintf("--------------------\n");
fprintf("Odd numbers from q_:\n");
odds_q_trip_3 = q_trip_3(q_trip_3_ind);
fprintf("%.2f\n", odds_q_trip_3);
fprintf("(odd - 1)/2:\n");
odds_q_trip_3 = double(odds_q_trip_3);
q_trip_3_num = (odds_q_trip_3 - ones(1, length(odds_q_trip_3)))/2;
fprintf("%.2f\n", q_trip_3_num);
fprintf("Number of triplets:\n");
fprintf("%.2f\n", sum(q_trip_3_num));



function IFFT_ = RIFFT(X)
    X = conj(X);
    n = length(X);
    x = RFFT(X)/n;
    IFFT_ = conj(x); % Divide by n to incorporate the normalization factor
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


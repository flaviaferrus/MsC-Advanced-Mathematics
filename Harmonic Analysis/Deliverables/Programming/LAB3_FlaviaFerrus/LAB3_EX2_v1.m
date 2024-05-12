%% LAB 3: FFT

%% EXERCISE 2

%N = 100000;
N = 7;
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
%fprintf("%.2f\n", suma);
fprintf("%.2f\n", p);


%p_fft = fft(p);
%q_fft = p_fft .* p_fft;
%q = ifft(q_fft);
%fprintf("q using FFT\n");
%fprintf("%.2f\n", q);
%fprintf("FFT of p:\n")
%fprintf("%f%+fj\n", real(p_fft), imag(p_fft));


%% Using regular polynomial multiplication:
q_2 = conv(p, p);
q_trip_3 = q_2(mod(q_2, 2) ~= 0);
q_trip_3_ind = q_trip_3 >= 3;
fprintf("--------------------\n");
fprintf("Odd numbers from q_:\n");
odds_q_trip_3 = q_trip_3(q_trip_3_ind);
fprintf("%.2f\n", odds_q_trip_3);
fprintf("(odd - 1)/2:\n");
q_trip_3_num = (odds_q_trip_3 - ones(1, length(odds_q_trip_3)))/2;
fprintf("%.2f\n", q_trip_3_num);
fprintf("Number of triplets:\n");
fprintf("%.2f\n", sum(q_trip_3_num));

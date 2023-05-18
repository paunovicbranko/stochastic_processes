close all
clear
clc
%% 

L = 1000;
N = 100;
p10 = 0.15;
p11 = 0.85;
p01 = 1 - p11;
p00 = 1 - p10;

X = zeros(L,N);
procena01 = zeros(1,L);

for i = 1:L
    %generisanje prvog clana pojedinacne realizacije
    X(i,1) = randi([0 1],1,1);
    
    for j = 2:N
        temp = rand(1,1);
        if X(i,j-1) == 1
            if temp < p11
                X(i,j) = 1;
            else
                X(i,j) = 0;
            end
        else
            if temp < p10
                X(i,j) = 1;
            else
                X(i,j) = 0;
            end
        end
    end
    procena01(i) = sum(X(i,1:end-1) == 0 & X(i,2:end) == 1) / (N-1);
end

teorijska_srednja_vrednost = 0.5 * p10;
srednja_vrednost_procena = sum(procena01) / length(procena01);
razlika_srednjih_vrednosti = abs(teorijska_srednja_vrednost - srednja_vrednost_procena);

varijansa_procena = sum((procena01-srednja_vrednost_procena).^2) / (length(procena01)-1);

[n, edges] = histcounts(procena01);
bin_width = edges(2) - edges(1);
pdf = n / sum(n);
figure,bar(edges(1:end-1), pdf, 'histc');
title('Normalizovani histogram');
xlabel('X osa');
ylabel('Verovatnoća');

%% 

%____________________________________________________________________________
%b
P10 = 0.85;
P11 = 0.85;
P00 = 1 - P10;
P01 = 1 - P11;
N1 = 100000;
new_array = zeros(1,N1);
new_array(1) = randi([0 1],1,1);
for j = 2:N1
    temp = rand(1,1);
    if new_array(j-1) == 1
        if temp < P11
            new_array(j) = 1;
        else
            new_array(j) = 0;
        end
    else
        if temp < P10
            new_array(j) = 1;
        else
            new_array(j) = 0;
        end
    end
end
m_b = 2 * N1 - 1;
rxx = zeros(1,m_b);
for i = 1:N1
    rxx(i) = new_array(N1-i+1:N1) * new_array(1:i)';
    rxx(m_b-i+1) = rxx(i)';
end
xosa_b = -25:25;
figure,stem(xosa_b,rxx(-25+N1:25+N1)/sum(new_array));
figure,stem(new_array(1:100));

%% 
%____________________________________________________________________
%c

Hmax = 0.5 * log2(1/0.5) + 0.5 * log2(1/0.5);
H0 = sum(new_array)/N1 * log2(1/(sum(new_array)/N1)) + (N1-sum(new_array))/N1 * log2(1/((N1-sum(new_array))/N1));

temporary_matrix = zeros(2);
for i = 2 : N1
    temporary_matrix(new_array(i-1) + 1, new_array(i)+1) = temporary_matrix(new_array(i-1)+1, new_array(i)+1)+1;
end

pr11 = temporary_matrix(1,1) / sum(temporary_matrix(1,:));
pr10 = temporary_matrix(1,2) / sum(temporary_matrix(1,:));
pr01 = temporary_matrix(2,1) / sum(temporary_matrix(2,:));
pr00 = temporary_matrix(2,2) / sum(temporary_matrix(2,:));

H1 = pr01 * log2(1/pr01) + pr11 * log2(1/pr11) + pr10 * log2(1/pr10) + pr00 * log2(1/pr00);


%% 
%____________________________________________________________________
%2

N2 = 500000;
a=-1; b=1;
X1 = a + (b-a) * rand(1,N2);

[n_2, edges_2] = histcounts(X1);
bin_width_2 = edges_2(2) - edges_2(1);
pdf_2 = n_2 / sum(n_2);
figure,bar(edges_2(1:end-1), pdf_2, 'histc');
title('Normalizovani histogram');
xlabel('X osa');
ylabel('Verovatnoća');

X12 = zeros(1,N2);
m=20;
for i = 1:m
    X12 = X12 + (a + (b-a) * rand(1,N2));
end

sr_vr = sum(X12)/length(X12);
var_vr = sum((X12-sr_vr).^2) / (length(X12)-1);

xosa_2 = -m:0.01:m;
gauss = exp(-(xosa_2 - sr_vr).^2/(2*var_vr))/(sqrt(var_vr*2*pi));

[n_2b, edges_2b] = histcounts(X12);
bin_width_2b = edges_2b(2) - edges_2b(1);
pdf_2b = n_2b / (sum(n_2b)*bin_width_2b);
figure,bar(edges_2b(1:end-1), pdf_2b, 'histc');
title('Normalizovani histogram');
xlabel('X osa');
ylabel('Verovatnoća');
hold on
plot(xosa_2,gauss,'r','LineWidth',2);



%% 
%____________________________________________________________________

lambda = 2;
N3 = 100000;
exponential_array = zeros(1,N3);

for i = 1:100
    uniform_array = rand(1,N3);
    %inverzna funkcija funkcije eksponencijalne raspodele:
    exponential_array = exponential_array + (-log(1-uniform_array)/lambda);
    
end

[n_3, edges_3] = histcounts(exponential_array);
bin_width_3 = edges_3(2) - edges_3(1);
pdf_3 = n_3 / sum(n_3);
figure,bar(edges_3(1:end-1), pdf_3, 'histc');
title('Normalizovani histogram');
xlabel('X osa');
ylabel('Verovatnoća');
obj = fitdist(exponential_array' ,'Gamma');
plot(obj);

%% 
%___________________________________________________________________
%c

N_3a = 100000;
c_array_x = randn(1,N_3a);
c_array_y = randn(1,N_3a);
r = sqrt(c_array_y .^2 + c_array_x .^2);

[n_3b, edges_3b] = histcounts(r);
bin_width_3b = edges_3b(2) - edges_3b(1);
pdf_3b = n_3b / (sum(n_3b)*bin_width_3b);
figure,bar(edges_3b(1:end-1), pdf_3b, 'histc');
title('Normalizovani histogram');
xlabel('X osa');
ylabel('Verovatnoća');
hold on

xosa_3b = 0:0.1:5;
rayleigh = xosa_3b .* exp(-xosa_3b.^2/2);

plot(xosa_3b, rayleigh,'LineWidth',2);

%% 

%___________________________________________________________________
%3

file_fausto = fopen('Fausto.txt','r');
data = fread(file_fausto);
figure,histogram(data,'Normalization','probability','BinWidth',1);

mean_value = mean(data);
var_value = var(data);
mean_square_value = meansqr(data);
meadiana = median(data);
modus = mode(data);

%autokovarijansa
[autocov, latency] = xcov(data, 20, 'normalized');
figure,stem(latency, autocov);
%___________________________________________________________________
%4
%% 

length_4 = 500000;
prob10 = 0.85;
prob11 = 0.85;
prob01 = 1 - prob11;
prob00 = 1 - prob10;

array_4 = zeros(1,length_4); %polazni niz

%generisanje prvog clana pojedinacne realizacije
array_4(1) = randi([0 1],1,1);

for j = 2:length_4
    temp = rand(1,1);
    if array_4(j-1) == 1
        if temp < prob11
            array_4(j) = 1;
        else
            array_4(j) = 0;
        end
    else
        if temp < prob10
            array_4(j) = 1;
        else
            array_4(j) = 0;
        end
    end
end

help_array = zeros(1,5*length_4); %pocetni niz gde je svaki clan ponovljen 5 puta

j = 1;
for i = 1:length_4
    if array_4(i) == 0
        help_array(j:j+4) = -1;
        j = j+5;
    else
        help_array(j:j+4) = 1;
        j = j+5;
    end
end

prob_of_error = zeros(1,11);

for k = 0 : 10
    
    snr = 10^(k/10);
    sigma = sqrt(1/snr);
    i_i_d_array = sqrt(sigma^2) * randn(1,5*length_4);
    reciever_input_array = i_i_d_array + help_array;
    
    %P(0) = P(0/0) * P(0) + P(0/1) * P(1)
    %P(0) = P(0/0) * P(0) + (1 - P(1/0)) * (1 - P(0))
    %P(0) = (1 - P(1/0)) / (2 - P(0/0) - P(1/0))
    
    prob0 = (1 - prob10) / (2 - prob00 - prob10);
    bopt = sigma^2 / 2 * log(prob0/(1-prob0));
    
    reciever_output_array = zeros(1,5 * length_4);
    
    for i = 1 : 5*length_4
        if reciever_input_array(i) <= bopt
            reciever_output_array(i) = -1;
        else
            reciever_output_array(i) = 1;
        end
    end
    
    j = 1;
    help_array_output = zeros (1, length_4);

    for i = 1 : 5 : 5*length_4
        if sum(reciever_output_array(i:i+4)) >= 0
            help_array_output(j) = 1;
        else
            help_array_output(j) = 0;
        end
        j = j + 1;
    end
   
    %prebrojava koliko ima razlika izmedju ulaza i izlaza
    bit_error_rate = sum(array_4 ~= help_array_output);
    prob_of_error(k+1) =bit_error_rate/length_4 ;

end

figure,semilogy(0:10,prob_of_error,'-o', 'LineWidth',2); hold on

p_error_map= 0.5 * erfc(sqrt((10.^((0:10)/10))/2)); 
neka_formula = nchoosek(5,3).*p_error_map.^3.*(1-p_error_map).^2+nchoosek(5,4).*p_error_map.^4.*(1-p_error_map)+nchoosek(5,5).*p_error_map.^5;

semilogy(0:10, neka_formula,'LineWidth',2);
title("Zavisnost verovatnoće greške od SNR-a");
xlabel('SNR [dB]');
ylabel('Verovatnoca greške [log]')
legend('Simulacija','Teorija');
%%
%_________________________________________________________
%b

prob_of_error_b = zeros(1,11);

for k = 0:10

    snr_b = 10^(k/10);
    sigma_b = sqrt(1/snr_b);
    i_i_d_array_b = sqrt(sigma_b^2) * randn(1,5*length_4);
    reciever_input_array_b = i_i_d_array_b + help_array;
    
    help_array_output_b = zeros (1, length_4);
    j = 1;
    
    for i = 1 : 5 : 5*length_4
        help_array_output_b(j) = sum(reciever_input_array_b(i:i+4));
        j = j+1;
    end
    
    reciever_output_array_b = zeros(1,length_4);
    
    for i = 1:length_4
        if help_array_output_b(i) <= bopt
            reciever_output_array_b(i) = 0;
        else
            reciever_output_array_b(i) = 1;
        end
    end
    
    bit_error_rate_b = sum(reciever_output_array_b ~= array_4);
    prob_of_error_b(k+1) = bit_error_rate_b/length_4;

end
p_error_map_b = 0.5*erfc(sqrt(5*(10.^((0:10)/10))/2));
figure,semilogy(0:10,prob_of_error_b,'-o','LineWidth',2);hold on
semilogy(0:10,p_error_map_b, 'LineWidth',2);
title("Zavisnost verovatnoće greške od SNR-a");
xlabel('SNR [dB]');
ylabel('Verovatnoca greške [log]')
legend('Simulacija','Teorija');
%%
%_________________________________________________________________________
%5

length_5 = 50000;
prob10_5 = 0.15;
prob11_5 = 0.85;
prob01_5 = 1 - prob11_5;
prob00_5 = 1 - prob10_5;

array_5 = zeros(1,length_5); %polazni niz

%generisanje prvog clana pojedinacne realizacije
array_5(1) = randi([0 1],1,1);

for j = 2:length_5
    temp = rand(1,1);
    if array_5(j-1) == 1
        if temp < prob11_5
            array_5(j) = 1;
        else
            array_5(j) = 0;
        end
    else
        if temp < prob10_5
            array_5(j) = 1;
        else
            array_5(j) = 0;
        end
    end
end

constelation = pskmod(0:3,4,pi/4,'gray');
binary_couple = reshape(array_5,2,length_5/2)';
symbols = zeros(length_5/2,1);

for i=1:length_5/2
    if binary_couple(i,1) == 0 && binary_couple(i,2) == 0
        symbols(i) = 0;
    else if  binary_couple(i,1) == 0 && binary_couple(i,2) == 1
            symbols(i) = 1;
    else if binary_couple(i,1) == 1 && binary_couple(i,2) == 0
            symbols(i) = 2;
    else 
        symbols(i) = 3;
    end
    end
    end
end

Signal(1:length_5/2,1) = constelation(symbols+1);

prob_of_error_5 = zeros(1,11);

for k = 0 : 10
    snr_5 = 10^((k+3)/10);
    sigma_5 = sqrt(1/snr_5);
    noise = complex(sigma_5 * randn(1,length_5/2),sigma_5 * randn(1,length_5/2));
    signal_reciever_input = Signal + noise.';

    if k == 1 || k == 7
        faze = [pi/4,3*pi/4,5*pi/4,7*pi/4]';
        figure,plot(real(signal_reciever_input(1:5000)),imag(signal_reciever_input(1:5000)),'ro');hold on
        scatter(cos(faze),sin(faze),'blue');
    end
        
    real_part = real(signal_reciever_input);
    imag_part = imag(signal_reciever_input);
    
    demodulated_signal = zeros(1,length_5);
    j=1;
    for i=1:length_5/2
        if real_part(i)<0 && imag_part(i)<0
            demodulated_signal(j) = 1;
            demodulated_signal(j+1) = 1;
            j=j+2;
        else if real_part(i)>0 && imag_part(i)<0
            demodulated_signal(j) = 1;
            demodulated_signal(j+1) = 0;
            j=j+2;
        else if real_part(i)<0 && imag_part(i)>0
            demodulated_signal(j) = 0;
            demodulated_signal(j+1) = 1;
            j=j+2;
        else 
            demodulated_signal(j) = 0;
            demodulated_signal(j+1) = 0;
            j=j+2;
        end
        end
        end
    end
    
    bit_error_rate_5 = sum(demodulated_signal~=array_5);
    prob_of_error_5(k+1) = bit_error_rate_5/length_5;
end

p_error_map_5 = 0.5*erfc(sqrt((10.^((0:10)/10))/2));
figure,semilogy(0:10,prob_of_error_5,'-o','LineWidth',2);hold on
semilogy(0:10 ,p_error_map_5,'LineWidth',2);
title("Zavisnost verovatnoće greške od SNR-a");
xlabel('SNR [dB]');
ylabel('Verovatnoca greške [log]')
legend('Simulacija','Teorija');



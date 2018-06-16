%ECE 257B
%Nonlinear parameter estimator
%NOTE: will take ~3minutes to run
%Paul Yoon
clear all; close all;

%CONTROLS (Yes=1, No=0)
ApplyNLAmplifier = 1;           %Use Nonlinear Amplifier?
ApplyWGN = 1;                   %Add white Gaussian noise?
ApplyMultipath = 0;             %Add Multipath Channel?
ApplyChannelEq = 0;             %Apply Receiver Channel Equalization?
ApplyPhaseEq = 0;               %Apply Phase Correction?
ApplyNLCorrection = 1;          %Apply Non-linear Correction?
ApplyNLEstimation =1;           %Estimate a1, a3?

PLOTCONSTELLATION = 0;          %Plot Constellations?
%Channel parameters
Gain=10;                        %Amplifier gain (dB)
IP3=-3;                         %Amplifier Nonlinearity parameter(dBm) 
snrDB_vec = 15:1:35;            %Signal-to-Noise Ratio (dB)
h = [1 0 0 -0.5];               %Channel Impulse Response (multipath channel)

ber_NLcomp = zeros(1,length(snrDB_vec));
ber_NL = zeros(1,length(snrDB_vec));



%Waveform params
MOD_ORDER               =  4;  
N_OFDM_SYMS             = 500;         % Number of OFDM symbols
SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
pilots = [1 1 -1 1].';
N_SC                    = 64;                                     % Number of subcarriers
CP_LEN                  = 16;                                     % Cyclic prefix length
N_DATA_SYMS             = N_OFDM_SYMS * length(SC_IND_DATA);      % Number of data symbols


%% TRANSMITTER
%preamble
   % LTS
   lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
   lts_t = ifft(lts_f, 64);
%    preamble = [lts_t lts_t zeros(1,32)]; %only lts preamble
LTS_reps = 20;
   preamble = [repmat([lts_t lts_t], 1, LTS_reps) zeros(1,32)]; %10 LTS
   for loops = 1:20
for m = 1:length(snrDB_vec)
snrdB = snrDB_vec(m)          %calls out the SNR vector at each element
%data
   number_of_bits= (N_DATA_SYMS * MOD_ORDER);
   tx_data = randi(2, 1, number_of_bits) - 1; 
   tx_syms = mapping(tx_data', MOD_ORDER, 1);
   tx_syms_mat = reshape(tx_syms, length(SC_IND_DATA), N_OFDM_SYMS);
%Pilots
   pilots_mat = repmat(pilots, 1, N_OFDM_SYMS);
%IFFT
   ifft_in_mat = zeros(N_SC, N_OFDM_SYMS);
   ifft_in_mat(SC_IND_DATA, :)   = tx_syms_mat; %comment out to null data
   ifft_in_mat(SC_IND_PILOTS, :) = pilots_mat;
   tx_payload_mat = ifft(ifft_in_mat, N_SC, 1);
   if(CP_LEN > 0)
      tx_cp = tx_payload_mat((end-CP_LEN+1 : end), :);
      tx_payload_mat = [tx_cp; tx_payload_mat];
   end
   tx_payload_vec = reshape(tx_payload_mat, 1, numel(tx_payload_mat));
   tx_vec = [preamble tx_payload_vec];

%% CHANNEL
%nonlinear amplifier
%corresponding taylor coefficients
a1 = 10^(Gain/10);
a3 = 2/(3*50)*10^((Gain-IP3)/10+3)/50; %50ohm system
tx_vec_L = a1*tx_vec;
if ApplyNLAmplifier == 1
    tx_vec_NL = a1*tx_vec + a3*tx_vec.*tx_vec.^2;

else
    tx_vec_NL = tx_vec_L;
end
%AWGN
snr = 10^(snrdB/10);
noisepower = var(tx_vec_NL)/(snr);
noise_vec = sqrt(noisepower)*(randn(1,length(tx_vec_NL)) + 1i*randn(1,length(tx_vec_NL)));
if ApplyWGN == 1
    tx_vec_NL = tx_vec_NL+noise_vec;
end
%channel impulse response
if ApplyMultipath == 1
    tx_vec_NL = filter(h,1,tx_vec_NL);
end

%% RECEIVER

%NL estimation (Least Squares Estimator)
rx_lts_avg = zeros(LTS_reps,64);
for i=0:(LTS_reps-1) %average preambles to reduce noise for NL estimate
    rx_lts_avg(i+1,:) = tx_vec_NL(i*length(lts_t)+1 : i*length(lts_t)+length(lts_t));
end
rx_lts_t_avg = sum(rx_lts_avg)/LTS_reps;
%{
    rx_lts_t_1 = tx_vec_NL(1:length(lts_t));  %extract received preamble observation vectors
    rx_lts_t_2 = tx_vec_NL(length(lts_t)+1:2*length(lts_t));
    rx_lts_t_avg = (rx_lts_t_1+rx_lts_t_2)/2;
%}
if ApplyNLEstimation ==1
    X = [lts_t.' (lts_t.^3).'];                 %Vandermonde Matrix
    a_est = inv(X'*X)*X'*rx_lts_t_avg.';         %least squares estimator
    a_est = real(a_est);                      %remove complex part
    a1_est = a_est(1);
    a3_est = a_est(2);
else
    a1_est = a1;
    a3_est = a3;
end
a1_est_vec(loops, m) = a1_est;
a3_est_vec(loops, m) = a3_est;
    tx_vec_NL = tx_vec_NL/a1;%normalized
% pkt detect and extract
    rx_preamble = tx_vec_NL(1:128);
    rx_payload_vec = tx_vec_NL(length(preamble)+1:end); %ideal case
% CFO est/correct
% remove cp
    rx_data_CFO = tx_vec_NL(length(preamble)+1:end);
    rx_payload_mat = reshape(rx_data_CFO,N_SC+CP_LEN,N_OFDM_SYMS);
    rx_payload_mat = rx_payload_mat(CP_LEN+1:end,:);
% FFT
    rx_payload_mat_fft = fft(rx_payload_mat,N_SC,1);
    
% Channel est/eq
%create nonlinear preamble to use in channel estimation so that channel
%equalization does not interfere with the nonlinear compensation
lts_NL = preamble + a3_est/a1_est*preamble.*preamble.^2;
lts_t_NL = lts_NL(65:128);
lts_f_NL = fft(lts_t_NL, N_SC);


rx_lts_t_cfo1 = rx_preamble(1:64);      %extract training sequences from preamble
rx_lts_t_cfo2 = rx_preamble(64+1:128);

rx_lts_f_cfo1 = fft(rx_lts_t_cfo1,64);              %take FFT of training sequences
rx_lts_f_cfo2 = fft(rx_lts_t_cfo2,64);
rx_lts_f_avg = (rx_lts_f_cfo1+rx_lts_f_cfo2)/2;
H_est = rx_lts_f_avg./lts_f_NL;                        %Channel estimate 

G = 1./H_est.';                                     %equalizer
G(isinf(G))=1;                                      %if a value of H is 0, then G is inf, instead use 1
% plot(abs(G))%testing
G_mat = repmat(G,1,N_OFDM_SYMS);
if ApplyChannelEq == 1
    rx_payload_mat_fft_eq = rx_payload_mat_fft.*G_mat;
else
    rx_payload_mat_fft_eq = rx_payload_mat_fft;
end
% phase eq
rx_pilots = rx_payload_mat_fft_eq(SC_IND_PILOTS,:);
rx_pilots_temp = rx_pilots.*pilots_mat;
phase_error = angle(sum(rx_pilots_temp)); %sum all the pilot subcarriers in ofdm block, and take the angle to find the avg phase error
phase_error_mat = repmat(phase_error,N_SC,1); %size the phase error matrix to the size of the symbol matrix
if ApplyPhaseEq==1
    rx_payload_mat_fft_eq = rx_payload_mat_fft_eq.*exp(-1i*phase_error_mat); %apply phase error equalization
end
%constellation plot
    rx_syms_mat = rx_payload_mat_fft_eq(SC_IND_DATA,:);
    rx_syms = reshape(rx_syms_mat,1,numel(rx_syms_mat));
    if PLOTCONSTELLATION ==1
    scatterplot(rx_syms(:))
    title(['Constellation with Nonlinear Distortion (SNR = ' num2str(snrdB) 'dB)'])
    end
%Demodulation
rx_data_NL = demapper(rx_syms,MOD_ORDER,1); 



%% NONLINEAR COMPENSATION
if ApplyNLCorrection ==1
%map rx bits into symbols
X_est = mapping(rx_data_NL, MOD_ORDER,1);
%form matrix for ifft
X_est_mat = reshape(X_est, length(SC_IND_DATA), N_OFDM_SYMS);
X_est_matp = zeros(N_SC, N_OFDM_SYMS);
X_est_matp(SC_IND_DATA, :)   = X_est_mat;
X_est_matp(SC_IND_PILOTS, :) = pilots_mat;
%take IFFT to get time domain symbols
x_est_time_mat = ifft(X_est_matp,N_SC,1);
%add CP
x_cp = x_est_time_mat((end-CP_LEN+1 : end), :);
x_est_time_mat = [x_cp; x_est_time_mat];
%P/S
x_est_time_vec = reshape(x_est_time_mat,1,numel(x_est_time_mat));
%calculate amplifer noise factor (REQUIRES KNOWING AMPLIFIER PARAMETERS)
NL_noise = a3_est/a1_est*(x_est_time_vec).*(x_est_time_vec).^2;
%S/P
NL_noise_mat = reshape(NL_noise, N_SC+CP_LEN, N_OFDM_SYMS);
%remove CP
NL_noise_mat = NL_noise_mat(CP_LEN+1:end,:);
%FFT to get back to freq domain
NL_noise_mat_freq = fft(NL_noise_mat,N_SC,1);
%subtract distortion term off originally decoded symbols
NL_corrected = rx_payload_mat_fft_eq-NL_noise_mat_freq;
%constellation
NL_corr_syms_mat = NL_corrected(SC_IND_DATA,:);
NL_corr_syms = reshape(NL_corr_syms_mat,1,numel(NL_corr_syms_mat));
if PLOTCONSTELLATION == 1
scatterplot(NL_corr_syms(:))
title(['Constellation with Nonlinear Distortion Correction SNR = ' num2str(snrdB) 'dB'])
end
end
%% BER Calculation

if ApplyNLCorrection ==1
rx_data_NLcomp = demapper(NL_corr_syms,MOD_ORDER,1); 
[number_NLcomp(m),ber_NLcomp(loops,m)] = biterr(tx_data,rx_data_NLcomp);
end
[number_NL(m),ber_NL(loops,m)] = biterr(tx_data,rx_data_NL);

end

   end
a1_est_vec_avg = sum(a1_est_vec)/loops;
a3_est_vec_avg = sum(a3_est_vec)/loops;
a1_est_err = abs(a1_est_vec_avg-a1)/a1 *100;
a3_est_err = abs(a3_est_vec_avg-a3)/a3 *100;
figure
plot(snrDB_vec, a1_est_err)
hold on 
plot(snrDB_vec, a3_est_err)
title('average NL estimation error')
xlabel('SNR'); ylabel('error (%)');
legend('a1 error','a3 error')
grid;


ber_NLcomp_avg = sum(ber_NLcomp)/loops;
ber_NL_avg= sum(ber_NL)/loops;
if ApplyNLCorrection ==1
    figure
    semilogy(snrDB_vec,ber_NL_avg,'r');
    title(['BER vs snrdB, IP3 = ' num2str(IP3)  'dBm'])
    xlabel('SNR (dB)')
    ylabel('Bit-Error-Rate (BER)')
    hold on
    semilogy(snrDB_vec,ber_NLcomp_avg);
    legend('Nonlinear','Nonlinear w/ compensation')
    hold off
    grid;
else
    figure
    semilogy(snrDB_vec,ber_NL,'r');
    title(['BER vs snrdB, IP3 = ' num2str(IP3)  'dBm'])
    xlabel('SNR (dB)')
    ylabel('Bit-Error-Rate (BER)')
    grid;   
end
 


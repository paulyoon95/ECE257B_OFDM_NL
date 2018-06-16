%ECE 257B
%Paul Yoon, Chak Chan
clear; close all;

% SOM parameters
mod_table = [1+1i 1-1i -1+1i -1-1i];
mod_table16 = [-3-3i -1-3i 1-3i 3-3i -3-1i -1-1i 1-1i 3-1i -3+1i -1+1i 1+1i 3+1i -3+3i -1+3i 1+3i 3+3i];

% CONTROLS (Yes=1, No=0)
ApplyNLAmplifier = 1;           %Use Nonlinear Amplifier?
ApplyWGN = 1;                   %Add white Gaussian noise?
ApplyMultipath = 0;             %Add Multipath Channel? <<Do not turn on>>
ApplyChannelEq = 0;             %Apply Receiver Channel Equalization? <<Do not turn on>>
ApplyNNOptimize = 0;            %Apply Neural Network Optimization <<Do not turn on>>
ApplyPhaseEq = 0;               %Apply Phase Correction? <<Do not turn on>>
ApplyNLCorrection = 1;          %Apply Non-linear Correction?
ApplySOM = 0;                   %Apply SOM Neural Network <<Do not turn on>>
ApplyIQPlot = 0;                %Apply IQ Plot during iterative cancellation process <<Can turn on>>

iterations = 10;                %Define number of iterations to do nonlinear estimation and cancellation

%% Channel parameters
Gain=0;                         %Amplifier gain (dB)
IP3=0;                         %Amplifier Nonlinearity parameter(dBm) 
Pin_dBm = -6;
snrDB_vec = 10:2:26;            %Signal-to-Noise Ratio (dB)
h = [1 0.9 0 -0.2];             %Channel Impulse Response (multipath channel)

ber_NLcomp = zeros(iterations,length(snrDB_vec));     % initialize to measure BER with nonlinear cancellation
ber_NL = zeros(1,length(snrDB_vec));                  % initialize to measure BER without nonlinear cancellation
symbol_error_rate = zeros(iterations+1,1);            % initialize to measure symbol_error_rate vs iterations


%% Transmitter, Channel and Receiver Setup
for m = 1:length(snrDB_vec)
    for j = 1:(ApplySOM+1)
        snrdB = snrDB_vec(m)          %calls out the SNR vector at each element

        %Waveform params
        MOD_ORDER               =   4;  
        N_OFDM_SYMS             = 500;                                    % Number of OFDM symbols
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
           preamble = [lts_t lts_t zeros(1,32)]; %only lts preamble
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
        Pin = 10^((Pin_dBm-30)/10);
        tx_vec = sqrt(Pin)*(tx_vec/(sqrt(sum(abs(lts_t).^2))))*length(lts_t);

        a1 = 10^((Pin_dBm - 30 + Gain)/10);
        a3 = 2/(3)*10^((Pin_dBm - 30 + Gain - IP3)/10); %50ohm system

        if ApplyNLAmplifier == 1
            tx_vec_NL = a1*tx_vec + a3*tx_vec.*abs(tx_vec).^2;
        else
            tx_vec_NL = tx_vec;
        end
        %AWGN << currently enabled >>
        snr = 10^(snrdB/10);
        noisepower = var(tx_vec_NL)/(snr);
        noise_vec = sqrt(noisepower)*(randn(1,length(tx_vec_NL)) + 1i*randn(1,length(tx_vec_NL)));
        if ApplyWGN == 1
            tx_vec_NL = tx_vec_NL+noise_vec;
        end
        %channel impulse response << currently disabled >>
        if ApplyMultipath == 1
            tx_vec_NL = filter(h,1,tx_vec_NL);
        end
        %% RECEIVER
        % pkt detect and extract
            rx_preamble = tx_vec_NL(1:128);
            rx_payload_vec = tx_vec_NL(length(preamble)+1:end); % detection algorithm not employed, assumed working

        % CFO est/correct
        % remove cp
            rx_data_CFO = tx_vec_NL(length(preamble)+1:end);
            rx_payload_mat = reshape(rx_data_CFO,N_SC+CP_LEN,N_OFDM_SYMS);
            rx_payload_mat = rx_payload_mat(CP_LEN+1:end,:);
        % FFT
            rx_payload_mat_fft = fft(rx_payload_mat,N_SC,1);

        % Channel est/eq
        rx_lts_t_cfo1 = rx_preamble(1:64);      %extract training sequences from preamble
        rx_lts_t_cfo2 = rx_preamble(64+1:128);

        rx_lts_f_cfo1 = fft(rx_lts_t_cfo1,64);              %take FFT of training sequences
        rx_lts_f_cfo2 = fft(rx_lts_t_cfo2,64);
        rx_lts_f_avg = (rx_lts_f_cfo1+rx_lts_f_cfo2)/2;
        H_est = rx_lts_f_avg./lts_f;                        %Channel estimate 

        G = 1./H_est.';                                     %equalizer
        G(isinf(G))=1;                                      %if a value of H is 0, then G is inf, instead use 1

        
        if ApplyChannelEq == 1
            G_mat = repmat(G,1,N_OFDM_SYMS);
            rx_payload_mat_fft_eq = rx_payload_mat_fft.*G_mat;
        else
            % << currently using this approach, no fading compensation yet >>
            H_est = sqrt(Pin)*(1/(sqrt(sum(abs(lts_t).^2))))*length(lts_t);
            rx_payload_mat_fft_eq = rx_payload_mat_fft/(a1*H_est);
        end

        if ApplyNNOptimize == 1
        % Neural Optimization
            nn_x_i = real(rx_lts_f_cfo1);
            nn_t_i = real(lts_f);
            net_i = feedforwardnet(10);
            net_i = train(net_i,nn_x_i,nn_t_i);

            nn_x_q = imag(rx_lts_f_cfo1);
            nn_t_q = imag(lts_f);
            net_q = feedforwardnet(10);
            net_q = train(net_q,nn_x_q,nn_t_q);
            for(z = 1:N_OFDM_SYMS)
                rx_payload_mat_fft_eq(z) = net_i(real(rx_payload_mat_fft_eq(z))) + net_q(imag(rx_payload_mat_fft_eq(z)));
            end
        end

        % phase eq
        rx_pilots = rx_payload_mat_fft_eq(SC_IND_PILOTS,:);
        rx_pilots_temp = rx_pilots.*pilots_mat;
        phase_error = angle(sum(rx_pilots_temp)); %sum all the pilot subcarriers in ofdm block, and take the angle to find the avg phase error
        phase_error_mat = repmat(phase_error,N_SC,1); %size the phase error matrix to the size of the symbol matrix
        if ApplyPhaseEq==1
            rx_payload_mat_fft_eq = rx_payload_mat_fft_eq.*exp(-1i*phase_error_mat); %apply phase error equalization
        end
        % constellation plot
            rx_syms_mat = rx_payload_mat_fft_eq(SC_IND_DATA,:);
            rx_syms = reshape(rx_syms_mat,1,numel(rx_syms_mat));
            if(ApplyIQPlot)
                scatterplot(rx_syms(:))
                title(['Constellation with Nonlinear Distortion (SNR = ' num2str(snrdB) 'dB)']);
            end
        % Demodulation
        rx_data_NL = demapper(rx_syms,MOD_ORDER,1);

        
        %% NONLINEAR COMPENSATION
        
        % map rx bits into symbols
        for n = 1:iterations
        if n == 1
            X_est = mapping(rx_data_NL, MOD_ORDER,1);
            X_est_temp = mapping(rx_data_NL, MOD_ORDER,1);
        else
            X_est = mapping(rx_data_NLcomp, MOD_ORDER,1);
        end
        symbol_error_rate(n,m) = size(find(abs((tx_syms_mat - reshape(X_est, length(SC_IND_DATA), N_OFDM_SYMS))) > 1*sqrt(2)),1);
        %form matrix for ifft
        X_est_mat = reshape(X_est, length(SC_IND_DATA), N_OFDM_SYMS);
        X_est_matp = zeros(N_SC, N_OFDM_SYMS);
        X_est_matp(SC_IND_DATA, :)   = X_est_mat;
        X_est_matp(SC_IND_PILOTS, :) = pilots_mat;
        % take IFFT to get time domain symbols
        x_est_time_mat = ifft(X_est_matp,N_SC,1);
        % add CP
        x_cp = x_est_time_mat((end-CP_LEN+1 : end), :);
        x_est_time_mat = [x_cp; x_est_time_mat];
        % P/S
        x_est_time_vec = reshape(x_est_time_mat,1,numel(x_est_time_mat));
        % calculate amplifer noise factor (REQUIRES KNOWING AMPLIFIER PARAMETERS)
        x_est_NL = a1*x_est_time_vec+a3*(x_est_time_vec).*(x_est_time_vec).^2;
        % NL_noise = x_est_NL-a1*x_est_time_vec; %noise is the a3 term
        NL_noise = (a3/a1)*H_est^2*(x_est_time_vec).*abs(x_est_time_vec).^2;
        % S/P
        NL_noise_mat = reshape(NL_noise, N_SC+CP_LEN, N_OFDM_SYMS);
        % remove CP
        NL_noise_mat = NL_noise_mat(CP_LEN+1:end,:);
        % FFT to get back to freq domain
        NL_noise_mat_freq = fft(NL_noise_mat,N_SC,1);
        % subtract distortion term off originally decoded symbols
        NL_corrected = rx_payload_mat_fft_eq-NL_noise_mat_freq;
        % constellation
        NL_corr_syms_mat = NL_corrected(SC_IND_DATA,:);
        NL_corr_syms = reshape(NL_corr_syms_mat,1,numel(NL_corr_syms_mat)); 
        
        if(ApplyIQPlot)
            scatterplot(NL_corr_syms(:))
            title(['Constellation with Nonlinear Distortion Correction (SNR = ' num2str(snrdB) 'dB)']);
        end
        
        %% BER Calculation
        rx_data_NLcomp = demapper(NL_corr_syms,MOD_ORDER,1); 
        [number_NL,ber_NL(1,m)] = biterr(tx_data,rx_data_NL);
        [number_NLcomp,ber_NLcomp(n,m)] = biterr(tx_data,rx_data_NLcomp);

        if(n==iterations)
            symbol_error_rate(n+1,m) = size(find(abs((tx_syms_mat - reshape(X_est, length(SC_IND_DATA), N_OFDM_SYMS))) > 1*sqrt(2)),1)/N_OFDM_SYMS;
        end
        
        if(n==1 && ApplySOM)
        %% SOM Training
        if j == 1
            b1 = 1;
            b2 = 0*b1;
            b3 = 0.007*b1
            b4 = 0*b1;
            b5 = 0*b1;
            b6 = 0*b1;
            NL_corr_syms_mod = b1*NL_corr_syms + b2*NL_corr_syms.^2 + b3*NL_corr_syms.^3 + b4*NL_corr_syms.^4 + b5*NL_corr_syms.^5 + b6*NL_corr_syms.^6;
            train_rx_iq = [real(NL_corr_syms_mod) ; imag(NL_corr_syms_mod)];
            net = selforgmap([4 4],20,3,'gridtop','linkdist');
            net = train(net,train_rx_iq);
            label_train_16QAM = [repmat([-3 -1 1 3],1,4); -3 -3 -3 -3 -1 -1 -1 -1 1 1 1 1 3 3 3 3];
            y = net(label_train_16QAM)
            classes = vec2ind(y);
        elseif(j==2 && ApplySOM)
            NL_corr_syms_mod = b1*NL_corr_syms + b2*NL_corr_syms.^2 + b3*NL_corr_syms.^3 + b4*NL_corr_syms.^4 + b5*NL_corr_syms.^5 + b6*NL_corr_syms.^6;
            rx_sym = zeros(1,length(NL_corr_syms_mod));
            y1 = net([real(NL_corr_syms_mod) ; imag(NL_corr_syms_mod)]);
            rx_id = vec2ind(y1);
            
            for i=1:(N_DATA_SYMS)
                match = find(classes == rx_id(i));
                rx_sym(i) = mod_table16(match);
            end
            
            rx_data_SOM = demapper(rx_sym,MOD_ORDER,1);
            
            % calculate SOM-Approach BER
            [number_SOM,ber_SOM(m)] = biterr(tx_data,rx_data_SOM);
        end
        end
        
        end
        
%         if(j == 1)
%             rx_symb_NLcomp = mapping2(rx_data_NLcomp,MOD_ORDER,1);
%             rx_symb_NLcomp_target = mapping2(tx_data', MOD_ORDER, 1);
%             save('rx_symb_NLcomp_train_set.mat','rx_symb_NLcomp');
%             save('rx_symb_NLcomp_train_target.mat','rx_symb_NLcomp_target');
%         elseif (j == 2)
%             rx_symb_NLcomp = mapping2(rx_data_NLcomp,MOD_ORDER,1);
%             rx_symb_NLcomp_target = mapping2(tx_data', MOD_ORDER, 1);
%             save('rx_symb_NLcomp_test_set.mat','rx_symb_NLcomp');
%             save('rx_symb_NLcomp_test_target.mat','rx_symb_NLcomp_target');
%         end

    end
end

figure(11)
plot((0:iterations),transpose([ber_NL(1,1); ber_NLcomp(:,6)]));
title(['BER-NL vs snrdB, ' 'IP3 = ' num2str(IP3)  'dBm'])
xlabel('# of Iterations')
ylabel('BER')
hold off


figure(50)
plot((0:iterations),transpose(symbol_error_rate(:,4)));         % Plot of Symbol Error Rate for 6th SNR vector elemnt 
title(['Symbol Error Rate vs Number of Iterations (IP3 = ' ,num2str(IP3), 'dBm; SNR = ', num2str(snrDB_vec(4)), 'dB)'])
xlabel('Number of Iterations')
ylabel('SER')

figure(51)
semilogy(snrDB_vec,ber_NL(1,:),'b-o');
hold on
grid on
semilogy(snrDB_vec,ber_NLcomp(2,:) , 'b--*');
title(['BER vs SNR (IP3 = ' ,num2str(IP3), 'dBm, Pin = ', num2str(Pin_dBm), 'dBm)'])
legend('Without NL-Compensation', 'With NL Compensation')
xlabel('SNR (dB)')
ylabel('BER')
hold off

ber_NL
ber_NLcomp
%% BER Simulation for 2x2 MIMO-OFDM with Rayleigh Fading and ZF Detection
clear all; close all; clc;

%% Simulation Parameters
N = 64;                 % Number of OFDM subcarriers
CP = 16;                % Cyclic prefix length (N/4)
N_sym = 100;            % Number of OFDM symbols per SNR
snr_db = 0:2:20;        % SNR range in dB
nTx = 2;                % Number of transmit antennas
nRx = 2;                % Number of receive antennas
mod_types = {'QPSK', '8PSK'};
L = 10;                 % Channel tap length (must be <= CP)

%% Initialize BER results
ber_results = zeros(length(snr_db), length(mod_types));

for mod_idx = 1:length(mod_types)
    modulation = mod_types{mod_idx};
    
    % Modulation parameters
    if strcmp(modulation, 'QPSK')
        M = 4;
        bits_per_sym = 2;
    else % 8PSK
        M = 8;
        bits_per_sym = 3;
    end
    
    % Calculate total symbols and bits needed
    sym_per_ofdm = N * nTx;          % 128 symbols per OFDM symbol
    bits_per_ofdm = sym_per_ofdm * bits_per_sym;
    
    for snr_idx = 1:length(snr_db)
        bit_errors = 0;
        total_bits_sent = 0;
        
        for sym_idx = 1:N_sym
            %% ============== Transmitter ==============
            % Generate random bits
            bits = randi([0 1], bits_per_ofdm, 1);
            
            % Reshape correctly for modulation
            bits_reshaped = reshape(bits, bits_per_sym, []);
            
            % Modulate
            mod_syms = pskmod(bits_reshaped, M, pi/M, 'gray', 'InputType', 'bit');
            
            % Reshape into N subcarriers × nTx antennas
            tx_freq = reshape(mod_syms, N, nTx);
            
            % OFDM Modulation (IFFT)
            tx_time = ifft(tx_freq, N); %%(problem here)
            
            % Add cyclic prefix
            tx_time_cp = [tx_time(end-CP+1:end, :); tx_time];
            
            %% ============== Channel ==============
            % Generate time-domain channel matrix
            H_time = (randn(L, nRx, nTx) + 1i*randn(L, nRx, nTx))/sqrt(2*L);
            
            % Compute frequency-domain channel
            H_freq = fft(H_time, N, 1);
            
            % Apply channel in time domain
            rx_time = zeros(N+CP, nRx);
            for rx = 1:nRx
                for tx = 1:nTx
                    % Convolve with channel impulse response
                    conv_result = conv(tx_time_cp(:, tx), squeeze(H_time(:, rx, tx)));
                    rx_time(:, rx) = rx_time(:, rx) + conv_result(1:N+CP);
                end
            end
            
            % Calculate noise power
            signal_power = mean(abs(tx_time_cp(:)).^2);
            noise_power = signal_power / (10^(snr_db(snr_idx)/10));
            
            % Add AWGN noise
            noise = sqrt(noise_power/2) * (randn(size(rx_time)) + 1i*randn(size(rx_time)));
            rx_time = rx_time + noise;
            
            %% ============== Receiver ==============
            % Remove cyclic prefix
            rx_time_no_cp = rx_time(CP+1:end, :);
            
            % OFDM Demodulation (FFT)
            rx_freq = fft(rx_time_no_cp, N);
            
            % ZF Detection for each subcarrier
            rx_syms = zeros(N, nTx);
            for k = 1:N
                H_k = squeeze(H_freq(k, :, :));
                % Zero-Forcing detection
                rx_syms(k, :) = (pinv(H_k) * rx_freq(k, :).').';
            end
            
            % Demodulation
            rx_bits = pskdemod(rx_syms(:), M, pi/M, 'gray', 'OutputType', 'bit');
            
            % Calculate errors
            bit_errors = bit_errors + sum(bits ~= rx_bits);
            total_bits_sent = total_bits_sent + length(bits);
        end
        
        % Store BER
        ber_results(snr_idx, mod_idx) = bit_errors/total_bits_sent;
        fprintf('%s at SNR %d dB: BER = %.4f\n', modulation, snr_db(snr_idx), ...
                ber_results(snr_idx, mod_idx));
    end
end

%% Plot Results
figure;
semilogy(snr_db, ber_results(:,1), 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(snr_db, ber_results(:,2), 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('2x2 MIMO-OFDM BER Performance (Rayleigh Fading, ZF Detection)');
legend('QPSK', '8PSK');
axis([min(snr_db) max(snr_db) 1e-5 1]);
set(gca, 'YScale', 'log');

% Parameters
nSubcarriers = 64;           % Number of subcarriers
nSymbols = 100;              % Number of OFDM symbols
cpLen = 16;                  % Cyclic prefix length
modOrder = 4;                % Modulation order (QPSK)

% Generate random data
dataBits = randi([0 1], nSubcarriers*log2(modOrder), nSymbols);

% QPSK Modulation
dataSymbols = qammod(dataBits, modOrder, 'InputType', 'bit', 'UnitAveragePower', true);

% Reshape into OFDM symbol matrix (each column is an OFDM symbol)
ofdmSymbols = reshape(dataSymbols, nSubcarriers, nSymbols);

% IFFT (OFDM modulation)
ofdmTime = ifft(ofdmSymbols, nSubcarriers);

% Add Cyclic Prefix
ofdmTimeCp = [ofdmTime(end-cpLen+1:end, :); ofdmTime];

% Channel (AWGN)
snr = 20; % Signal-to-Noise Ratio in dB
rxSignal = awgn(ofdmTimeCp, snr, 'measured');

% Receiver
% Remove Cyclic Prefix
rxSignalCpRemoved = rxSignal(cpLen+1:end, :);

% FFT (OFDM demodulation)
rxSymbols = fft(rxSignalCpRemoved, nSubcarriers);

% QPSK Demodulation
rxDataBits = qamdemod(rxSymbols, modOrder, 'OutputType', 'bit', 'UnitAveragePower', true);

% Calculate BER (Bit Error Rate)
[~, ber] = biterr(dataBits(:), rxDataBits(:));
disp(['Bit Error Rate (BER): ' num2str(ber)]);

% Plot constellation diagram
scatterplot(rxSymbols(:));
title('Received Constellation Diagram');

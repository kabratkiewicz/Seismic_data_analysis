%% Lissajous signal analysis
% Maciej J. Mendecki, K. Abratkiewicz, "Rammstein Concert-Induced 
% Seismicity Signal Filtering by Time-Frequency Analysis", IEEE Geoscience
% and Remote Sensing Letters, 2024
% Data: Maciej J. Mendecki, maciej.mendecki@us.edu.pl, Faculty of Natural
% Sciences, University of Silesia in Katowice, 41-200 Sosnowiec, Poland.
% software: Karol Abratkiewicz, karol.abratkiewicz@pw.edu.pl, Institute 
% of Electronic Systems, Warsaw University of Technology, 00-665 Warsaw, 
% Poland

close all
clear
clc

%% initialization
fontsize = 20;    % fontsize
img_max_size = 1; % 1 - fullscreen, otherwise - default size
Init_Env(fontsize, img_max_size)

%% EW component filtering
load("data\sunday\Rammstein_concert_sunday_EW.mat") % load component
signal_EW  = signal .* 7.94E-10; % scaling to V/s;
d          = 8; % decimation factor
signal_EW  = decimate(signal_EW,d,"fir"); % signal decimation decimation
fs         = fs/d; % sampling rate
N          = length(signal_EW);
NFFT       = 2048; % number of FFT points (both positive and negative frequencies)
L          = 150; % Gaussian widnow time spread
thr        = -50; % threshold for spectrogram
method     = 'FFT'; % spectrogram computation method (FFT or ptByPt)

%% signal structure initialization and definition
x.signal = signal_EW; % ssignal vector
x.N      = length(signal_EW); % signal length
x.fs     = fs; % sampling rate

%% spectrogram computation
S        = Gab_STFT(x, NFFT, L, method); % STFT
f_scale  = linspace(-x.fs/2, x.fs/2, NFFT); % frequency scale
t_scale  = linspace(0,x.N/x.fs,x.N)./60;    % time scale
Plot_Energy(S, thr, 0, t_scale, f_scale, fontsize, img_max_size); % display spectrogram
ylim([0 fs/2])
drawnow

%% STFT masking
mask = GenMask(N, NFFT); % mask generation
rec_EW = Gab_ISSTFT(S .* (mask), L, NFFT); %signal reconstruction from ISTFT
rec_EW = real(rec_EW); % RECONSTRUCTED (FILTERED) SIGNAL!! only the real part


%% NS component filtering
load("data\sunday\Rammstein_concert_sunday_NS.mat") % load component
signal_NS  = signal .* 7.94E-10; % scaling to V/s;
d          = 8; % decimation factor
signal_NS  = decimate(signal_NS,d,"fir"); % signal decimation decimation
fs         = fs/d; % sampling rate
N          = length(signal_NS);
NFFT       = 2048; % number of FFT points (both positive and negative frequencies)
L          = 150; % Gaussian widnow time spread
thr        = -50; % threshold for spectrogram
method     = 'FFT'; % spectrogram computation method (FFT or ptByPt)

%% signal structure initialization and definition
x.signal = signal_NS; % ssignal vector
x.N      = length(signal_NS); % signal length
x.fs     = fs; % sampling rate

%% spectrogram computation
S        = Gab_STFT(x, NFFT, L, method); % STFT
f_scale  = linspace(-x.fs/2, x.fs/2, NFFT); % frequency scale
t_scale  = linspace(0,x.N/x.fs,x.N)./60;    % time scale
Plot_Energy(S, thr, 0, t_scale, f_scale, fontsize, img_max_size); % display spectrogram
ylim([0 fs/2])
drawnow

%% STFT masking
mask   = GenMask(N, NFFT); % mask generation
rec_NS = Gab_ISSTFT(S .* (mask), L, NFFT); %signal reconstruction from ISTFT
rec_NS = real(rec_NS); % RECONSTRUCTED (FILTERED) SIGNAL!! only the real part


%% Lissajous analysis
M       = 4; % analyzed signal lenght
NFrames = floor(N/M); % number of M-second frames
frame   = M * fs; % frame length
for i = 732 %1 : NFrames % fixed to show only a selected frame
    originalSignalPartEW = signal_EW((i-1)*frame + 1 : i * frame); % original EW signal part
    originalSignalPartNS = signal_NS((i-1)*frame + 1 : i * frame); % original NS signal part

    filteredSignalPartEW = rec_EW((i-1)*frame + 1 : i * frame);    % filtered EW signal part
    filteredSignalPartNS = rec_NS((i-1)*frame + 1 : i * frame);    % filtered NS signal part
    
    % plot the Lissajous curves
    figure
    plot(originalSignalPartEW.*1e6, originalSignalPartNS.*1e6,'r','LineWidth',3)
    hold on
    plot(filteredSignalPartEW.*1e6,filteredSignalPartNS.*1e6,'k','LineWidth',3)
    grid on
    grid minor
    axis equal
    legend('Before filtering','After filtering')
    % saveas(gca, append('lissajous/fig', num2str(i),'.png'))
    ylabel('A [$\mu $m/s]')
    xlabel('A [$\mu $m/s]')
    drawnow
end

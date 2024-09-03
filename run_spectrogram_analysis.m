%% Time-frequency filtering of seismic signals
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

%% parameters
load("data\sunday\Rammstein_concert_sunday_EW.mat") % load component
signal  = signal .* 7.94E-10; % scaling to V/s;
d       = 8; % decimation factor
signal  = decimate(signal,d,"fir"); % signal decimation decimation
fs      = fs/d; % sampling rate
N       = length(signal);
NFFT    = 2048; % number of FFT points (both positive and negative frequencies)
L       = 150; % Gaussian widnow time spread
thr     = -50; % threshold for spectrogram
method  = 'FFT'; % spectrogram computation method (FFT or ptByPt)

%% signal structure initialization and definition
x.signal = signal; % ssignal vector
x.N      = length(signal); % signal length
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
rec = Gab_ISSTFT(S .* (mask), L, NFFT); %signal reconstruction from ISTFT
rec = real(rec); % RECONSTRUCTED (FILTERED) SIGNAL!! only the real part

%% mask display
figure; 
imagesc(t_scale,f_scale,mask); set(gca,'ydir','normal')
ylim([0 fs/2])
colormap(flipud(inferno(256)))
xlabel('t [min]','FontSize',fontsize, 'interpreter','latex')
ylabel('f [Hz]','FontSize',fontsize, 'interpreter','latex')
%%
x.signal = awgn(rec, 200,'measured'); % add a very small noise - otherwise the STFT is computed from zeros
%% reconstructed signal display
x.N = length(rec);
x.fs = fs;
S = Gab_STFT(x, NFFT, L,method);
f_scale = linspace(-x.fs/2, x.fs/2, NFFT) ;
t_scale = linspace(0,x.N/x.fs,x.N)./60;
Plot_Energy(S , thr, 0, t_scale, f_scale, fontsize, img_max_size);
ylim([0 fs/2])
drawnow



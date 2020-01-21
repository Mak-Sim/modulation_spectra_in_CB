function [features, cfg] = MS_wdftfb_analysis(x, fs)
%WDFTFB_MS_ANALYSIS Analysis of signal using warped DFT-modulated filter bank
%                in bark frequency scale. Bark-scaled frequency
%                coefficeints are extracted. Modulation features are
%                computed.

T   = 800;    % ms (supra-frame size)
load('WDFT_FB_M40_h200.mat');

Nfb       = params.M/2;
R1 = fs/params.fs;            % downsampling factor 1
R2 = params.Sk(params.M/2);   % downsampling factor 2

Nframe    = round(T/1000*params.fs/R2);
Nfshift   = round(Nframe*0.25);

fs3 = 512;       % Hz (modulation frequency)
D = (fs/R1/R2)/fs3;

% Window functions
wnd_mod = hamming(ceil(Nframe/D))';

%1. Downsampling by R_1 (in the scheme)
x = resample(x,params.fs,fs);

%2. Warped DFT-modulated FB
X = warped_dftfb_analysis(x, params);
X = X(1:Nfb,length(params.h)*3:end);    % removing transient process in bands
X_D = X(:,1:R2:end);  % downsampling by R_2 (can be done inside filter bank stucture)

FBE = abs(X_D).^2;

%% Second transform (to modulation frequencies)
fs2 = round(params.fs/R2);
fs_mod = 512;       % Hz
frc_mod_high = fs_mod/2;  % Hz
N_frc_mod = 257;
frc_mod = linspace(0,frc_mod_high,N_frc_mod);

K = floor((size(FBE, 2) - Nframe)/Nfshift);

% features = zeros(Nfb*N_frc_mod, K); % 1-D feature vectors
features = zeros(Nfb,N_frc_mod, K); % 2-D feature vectors


for N = 1:K % number of frames
    id = 1+(N-1)*Nfshift : Nframe+(N-1)*Nfshift;
    X_frc_mod = zeros(Nfb,N_frc_mod);
    
    for k=1:Nfb
        k_fb_band = FBE(k,id);
        k_fb_band = k_fb_band - mean(k_fb_band);
        FBE_k_ds = resample(k_fb_band,fs_mod,fs2);
        X_frc_mod(k,:) = abs(freqz(FBE_k_ds.*wnd_mod,1,frc_mod,fs_mod));        
    end
%     features(:, N) = 10*log10(X_frc_mod(:));    % 1-D feature vectors
    features(:,:, N) = 10*log10(X_frc_mod);    % 2-D feature vectors
end

%% Storing of some parameters
cfg.supra_frame_size_ms = T;
cfg.CB_num = Nfb;       
cfg.frc_mod = frc_mod;          % modulation frequencies
end


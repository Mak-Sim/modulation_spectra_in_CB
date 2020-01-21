function [features, cfg] = MS_fft_analysis(x, fs)

T   = 800;    % ms (supra-frame size)
frame_size_ms = 40; % ms
h_size_ms = 0.47;   % ms (hop size)

Nfb       = 20; % Number of critical bands
frame_size = round(frame_size_ms/1000*fs);
h_size   = round(h_size_ms/1000*fs);
N_bins = floor(frame_size/2) + 1;

fs2 = round(fs/h_size);
fs3 = 512;       % Hz (modulation frequency)
D = fs2/fs3;

% Supra-frames parameters
Nframes   = round((T-frame_size_ms)/h_size_ms) + 1;     % number of frames in supra-frame
Nfshift   = round(Nframes*0.25);                   % overlapping of supra-frames

% window functions
wnd_acoust = hamming(frame_size);
wnd_mod = hamming(ceil(Nframes/D))';


Nframes_total = floor((length(x) - frame_size)/h_size);

%% This section is needed for matching this analysis scheme with scheme based on warped filter bank
fs_fb = 12600;
alpha = (0.1957 - 1.048*((2/pi)*atan(0.07212*(fs_fb/1000))).^(1/2)); % Selection of warping coefficient (Bark scale)
fb_edge = pi/(2*Nfb):pi/Nfb:pi; % Edges of the filter bank
w_def_edge = -freq_warp(fb_edge,alpha);
h_edge = round(w_def_edge/(2*pi)*fs_fb);
f_edge = interp1(0:fs/frame_size:fs,1:frame_size+1,h_edge);

H = zeros(Nfb, N_bins);

% Rectangular windows (in MFCC -- triangle windows)
for i = 1:Nfb
    if (i==1)
        H(i, :) = trapmf(1:N_bins,[1 1 f_edge(i) ceil(f_edge(i))]) ;    
    else
        H(i, :) = trapmf(1:N_bins, [floor(f_edge(i-1)) f_edge(i-1) f_edge(i) ceil(f_edge(i))]);
    end
end

% 1st Transform
X_frc_time = zeros(Nfb,Nframes_total);
for N=1:Nframes_total
    x_frame = x(1+(N-1)*h_size:(N-1)*h_size+frame_size);
    x_frame = x_frame.*wnd_acoust;
    X = fft(x_frame); 
    X = X(1:N_bins);
            
    %% MFCC-like calculation
    Y = abs(X).^2;
    FBE = H*Y;  
    
    %% Complex-signal summation
%     FBE = abs(H*X).^2;
    
    X_frc_time(:,N) = FBE;
end

% Second transform
frc_mod_high = fs3/2;  % Hz
frc_mod = 0:1:frc_mod_high;
N_frc_mod = length(frc_mod);

K = floor((Nframes_total - Nframes)/Nfshift);
% features = zeros(Nfb*N_frc_mod, K);   % 1-D feature vector
features = zeros(Nfb, N_frc_mod, K);     % 2-D feature vector

for N = 1:K % number of frames
    id = 1+(N-1)*Nfshift : Nframes+(N-1)*Nfshift;
    X_frc_mod = zeros(Nfb,N_frc_mod);
    
    for k=1:Nfb
        k_fb_band = X_frc_time(k,id);
        k_fb_band = k_fb_band - mean(k_fb_band);
        FBE_k_ds = resample(k_fb_band,fs3,fs2);
        X_frc_mod(k,:) = abs(freqz(FBE_k_ds.*wnd_mod,1,frc_mod,fs3)); 
    end
    
%     features(:, N) = 10*log10(X_frc_mod(:));   % 1-D feature vector
    features(:,:, N) = 10*log10(X_frc_mod);   % 2-D feature vector
end

%% Storing of some parameters
cfg.supra_frame_size_ms = T;
cfg.CB_num = Nfb;   
cfg.frame_size = frame_size;    
cfg.frc_mod = frc_mod;          % modulation frequencies
end
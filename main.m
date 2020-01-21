% Calculation of modulation spectra in critical bands. Two different
% approaches are used for estimation energy is critical bands: 
%   1) DFT with channel merging;
%   2) applaying warped DFT-modulated filter bank

addpath('auxilary');

[x,fs] = audioread('ent_h_020.wav');
% [x,fs] = audioread('p020.wav');

[MS_in_CB, cfg] = MS_fft_analysis(x, fs);
[MSw_in_CB, ~]  = MS_wdftfb_analysis(x, fs);

max_E = max(max(MS_in_CB(:,:,1)));
max_Ew = max(max(MSw_in_CB(:,:,1)));
[F_mod,CB] = meshgrid(cfg.frc_mod, 1:cfg.CB_num);

%%
figure;
subplot(2,2,[1 2]);
N = cfg.supra_frame_size_ms/1000*fs;
plot((1:N)/fs,x(1:N));
xlabel('�����{\it t}, {\it �}');
ylabel('\textit{x}(\textit{t})','interp','latex');
title('������� ������');
    
subplot(2,2,3)
pcolor(F_mod,CB,MS_in_CB(:,:,1));
shading flat;
colormap('jet');
set(gca,'CLim',[max_E-60 max_E]);

xlabel('������������� �������, {\it ��}');
ylabel('����� ����������� ������');
title({'������������� ������', '(��� � ������������ �����)'});

subplot(2,2,4)
pcolor(F_mod,CB,MSw_in_CB(:,:,1));
shading flat;
colormap('jet');
set(gca,'CLim',[max_Ew-50 max_Ew]);

xlabel('������������� �������, {\it ��}');
ylabel('����� ����������� ������');
title({'������������� ������','(��������������� ���� ��������)'});


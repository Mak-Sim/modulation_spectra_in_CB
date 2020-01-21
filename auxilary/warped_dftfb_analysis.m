function [ X ] = warped_dftfb_analysis(x, params)
%warped_cmfb_analysis - ��������� ���-�������������� ���� �������� �������.
%   x    --  ������������� ������.
%  alpha -- �������� ����������� �������

% 1. ���������
M = params.M;       % ����� ������� ����� ��������
N = params.N;       % ������� �������-��������� (��)
m = params.m;       % ������� ���������� ��������� ��
alpha = params.alpha;
h = params.h;       % ������������ �������-���������
h = h/sum(h);       % ����������

Npt = length(x);
X = zeros(M, Npt);                   % ������� ��������� ��������
allpass_delay_chain = zeros(N,1);
  
% ���������� ����� �������� �������
polyphase_sum = sparse(kron(ones(1,m),eye(M)));
for n=1:Npt
    % 1. ������� ������� �������
    R1 = x(n);
    R2 = allpass_delay_chain(1);
    allpass_delay_chain(1) = R1;
    for i=2:N
        R3 = R2;
        R2 = allpass_delay_chain(i);
        R1 = (R2-R1)*alpha + R3;
        allpass_delay_chain(i) = R1; 
    end
    
    % 2. ��������� �� ������������ �������-���������
    PP_out = allpass_delay_chain.*h;
    
    Xk = polyphase_sum*PP_out;

    % 3. ���������� ��������� 
    X(:,n) = M*ifft(Xk);
end
end


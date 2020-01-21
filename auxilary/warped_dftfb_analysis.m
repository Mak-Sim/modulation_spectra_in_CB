function [ X ] = warped_dftfb_analysis(x, params)
%warped_cmfb_analysis - реализуют ДПФ-модулированный банк фильтров анализа.
%   x    --  анализируемый сигнал.
%  alpha -- параметр искривления частоты

% 1. Параметры
M = params.M;       % Число каналов банка фильтров
N = params.N;       % Порядок фильтра-прототипа (ФП)
m = params.m;       % Порядок полифазных компонент ФП
alpha = params.alpha;
h = params.h;       % Коэффициенты фильтра-прототипа
h = h/sum(h);       % Нормировка

Npt = length(x);
X = zeros(M, Npt);                   % матрица канальный сигналов
allpass_delay_chain = zeros(N,1);
  
% Реализация банка фильтров анализа
polyphase_sum = sparse(kron(ones(1,m),eye(M)));
for n=1:Npt
    % 1. Цепочка фазовых звеньев
    R1 = x(n);
    R2 = allpass_delay_chain(1);
    allpass_delay_chain(1) = R1;
    for i=2:N
        R3 = R2;
        R2 = allpass_delay_chain(i);
        R1 = (R2-R1)*alpha + R3;
        allpass_delay_chain(i) = R1; 
    end
    
    % 2. Умножение на коэффициенты фильтра-прототипа
    PP_out = allpass_delay_chain.*h;
    
    Xk = polyphase_sum*PP_out;

    % 3. Косинусная модуляция 
    X(:,n) = M*ifft(Xk);
end
end


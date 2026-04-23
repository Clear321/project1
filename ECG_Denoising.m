%% ECG信号降噪处理程序
% =========================================================================
% 基于MATLAB的数字滤波器设计，用于心电图（ECG）信号降噪处理。
% 针对ECG信号中常见的高斯白噪声、工频干扰及基线漂移三种噪声类型，
% 分别设计IIR与FIR数字滤波器，以有效抑制噪声干扰，提升信号质量。
%
% 滤波器设计指标（采样率 Fs = 500 Hz）：
%   高斯白噪声滤除（低通）：fp = 25 Hz,  fs = 40 Hz,  Ap = 1 dB, As = 15 dB
%   工频干扰滤除  （带阻）：fp = [45,55] Hz, fs = [48,52] Hz, Ap = 1 dB, As = 15 dB
%   基线漂移滤除  （高通）：fp = 0.7 Hz, fs = 0.3 Hz, Ap = 1 dB, As = 15 dB
% =========================================================================

clear; clc; close all;

%% =========================================================================
%  全局参数
%% =========================================================================
Fs = 500;   % 采样率 (Hz)
Ap = 1;     % 通带最大衰减 (dB)
As = 15;    % 阻带最小衰减 (dB)

%% =========================================================================
%  第0节：导入或生成ECG信号
%% =========================================================================
if exist('Ecg.txt', 'file')
    val = importdata('Ecg.txt');   % 从文件读取真实ECG数据
    % importdata 返回结构取决于文件格式：纯数值文件返回矩阵，带标题的文件返回struct
    if isstruct(val)
        raw = val.data;
    else
        raw = val;
    end
    % 兼容行向量和矩阵：优先取第一行，不足1800点则截断或补全
    raw = raw(:)';            % 统一转为行向量
    if length(raw) >= 1800
        signal = raw(1:1800);
    else
        signal = [raw, zeros(1, 1800 - length(raw))];
        warning('Ecg.txt数据不足1800点，已用零补全。');
    end
else
    % Ecg.txt不存在时，自动生成合成ECG信号（模拟真实心电图波形）
    signal = gen_ecg(1800, Fs);
end

N  = length(signal);            % 信号长度
t0 = (0:N-1) / Fs;             % 时间轴 (s)
f  = (0:N-1) * (Fs / N);       % 频率轴 (Hz)

%% =========================================================================
%  第1节：生成含三种噪声的ECG信号
%% =========================================================================

% 1.1 加入高斯白噪声（SNR = 10 dB）
%     awgn()按照信号实际功率加入高斯白噪声
Signal1 = awgn(signal, 10, 'measured');

% 1.2 加入50Hz工频干扰（余弦信号）
%     工频干扰主要由电力系统引起
av      = 100;                              % 干扰幅值
t_idx   = 1:N;
noise2  = av * cos(2*pi*50*t_idx/Fs);      % 50 Hz 余弦干扰
Signal2 = signal + noise2;

% 1.3 加入基线漂移噪声（低频线性漂移）
%     模拟呼吸、肢体运动等引起的基线漂移
n1      = floor(N/3);
x1      = zeros(1, n1);                     % 前1/3段无漂移
t_bl    = 1:(N - n1);
x2      = (N - n1) / 2000 * (t_bl - 1) + 1;  % 后2/3段线性漂移
noise3  = [x1, x2];
Signal3 = signal + noise3;

%% =========================================================================
%  Figure 1：四路信号时域波形对比
%% =========================================================================
figure(1);
subplot(4,1,1); plot(t0, signal, 'k');
title('干净的ECG信号'); xlabel('时间/s'); ylabel('幅值'); grid on;

subplot(4,1,2); plot(t0, Signal1, 'k');
title('加入高斯白噪声的信号 (SNR=10dB)'); xlabel('时间/s'); ylabel('幅值'); grid on;

subplot(4,1,3); plot(t0, Signal2, 'k');
title('加入工频干扰的信号 (50Hz)'); xlabel('时间/s'); ylabel('幅值'); grid on;

subplot(4,1,4); plot(t0, Signal3, 'k');
title('加入基线漂移的信号'); xlabel('时间/s'); ylabel('幅值'); grid on;

sgtitle('ECG信号与三种噪声 — 时域波形对比');

%% =========================================================================
%  Figure 2：四路信号频域幅度谱对比
%% =========================================================================
figure(2);
plot_four_spectra(signal, Signal1, Signal2, Signal3, N, f);
sgtitle('ECG信号与三种噪声 — 频域幅度谱对比');

%% =========================================================================
%  第2节：高斯白噪声滤除 — 低通滤波
%  设计指标：fp = 25 Hz，fs = 40 Hz，Ap = 1 dB，As = 15 dB，Fs = 500 Hz
%% =========================================================================
fp_lp = 25;    % 低通通带截止频率 (Hz)
fs_lp = 40;    % 低通阻带截止频率 (Hz)

% 数字归一化频率（相对于奈奎斯特频率 Fs/2）
Wp_lp = fp_lp / (Fs/2);
Ws_lp = fs_lp / (Fs/2);

%% ---- 2.1 低通IIR滤波器（Butterworth） ----
% buttord()：根据设计指标计算最小所需阶数和-3dB截止频率
% butter()：设计Butterworth低通数字滤波器
[ord_iir_lp, Wn_iir_lp] = buttord(Wp_lp, Ws_lp, Ap, As);
[b_iir_lp,   a_iir_lp ] = butter(ord_iir_lp, Wn_iir_lp, 'low');

% filter()：对含噪信号进行IIR滤波
Sig1_IIR = filter(b_iir_lp, a_iir_lp, Signal1);

% Figure 3：IIR低通滤波器幅频特性 + 滤波前后时域对比
figure(3);
subplot(2,1,1);
[h3, w3] = freqz(b_iir_lp, a_iir_lp, 2048, Fs);
plot(w3, 20*log10(abs(h3)+eps), 'k', 'LineWidth', 1.2);
title(sprintf('低通IIR（Butterworth）幅频特性  阶数 N = %d', ord_iir_lp));
xlabel('频率/Hz'); ylabel('幅值/dB'); xlim([0 Fs/2]); grid on;
xline(fp_lp, '--b', sprintf('f_p=%d Hz', fp_lp));
xline(fs_lp, '--r', sprintf('f_s=%d Hz', fs_lp));

subplot(2,1,2);
plot(t0, Signal1, 'Color', [0.6 0.6 0.6]); hold on;
plot(t0, Sig1_IIR, 'k', 'LineWidth', 1.2); hold off;
legend('含高斯白噪声', 'IIR低通滤波后', 'Location', 'best');
title('高斯白噪声 — IIR低通滤波效果（时域）');
xlabel('时间/s'); ylabel('幅值'); grid on;
sgtitle('高斯白噪声滤除 — 低通IIR（Butterworth）滤波器');

% Figure 4：IIR低通滤波前后频域对比
figure(4);
compare_spectra(Signal1, Sig1_IIR, N, f, Fs, ...
    '含高斯白噪声信号频谱', 'IIR低通滤波后频谱', ...
    '高斯白噪声滤除 — 低通IIR频域对比');

%% ---- 2.2 低通FIR滤波器（Hamming窗） ----
% 利用Hamming窗经验公式估算滤波器阶数：N ≈ 6.6 / (Δf/Fs)
% 其中 Δf 为过渡带宽（Hz）
delta_f_lp = (fs_lp - fp_lp) / Fs;     % 归一化过渡带宽
ord_fir_lp = ceil(6.6 / delta_f_lp);   % 估算FIR阶数

% fir1()要求高通/带阻时阶数为偶数（滤波器长度为奇数），低通两者皆可
% 为统一，此处低通也保持偶数阶
if mod(ord_fir_lp, 2) ~= 0
    ord_fir_lp = ord_fir_lp + 1;
end

% 截止频率取通带与阻带截止频率的中点
Wc_fir_lp = (fp_lp + fs_lp) / 2 / (Fs/2);
b_fir_lp  = fir1(ord_fir_lp, Wc_fir_lp, 'low', hamming(ord_fir_lp+1));

% fftfilt()：使用FFT重叠相加法对信号进行FIR滤波（适合高阶FIR）
Sig1_FIR = fftfilt(b_fir_lp, Signal1);

% Figure 5：FIR低通滤波器幅频特性 + 滤波前后时域对比
figure(5);
subplot(2,1,1);
[h5, w5] = freqz(b_fir_lp, 1, 2048, Fs);
plot(w5, 20*log10(abs(h5)+eps), 'k', 'LineWidth', 1.2);
title(sprintf('低通FIR（Hamming窗）幅频特性  阶数 N = %d', ord_fir_lp));
xlabel('频率/Hz'); ylabel('幅值/dB'); xlim([0 Fs/2]); grid on;
xline(fp_lp, '--b', sprintf('f_p=%d Hz', fp_lp));
xline(fs_lp, '--r', sprintf('f_s=%d Hz', fs_lp));

subplot(2,1,2);
plot(t0, Signal1, 'Color', [0.6 0.6 0.6]); hold on;
plot(t0, Sig1_FIR, 'k', 'LineWidth', 1.2); hold off;
legend('含高斯白噪声', 'FIR低通滤波后', 'Location', 'best');
title('高斯白噪声 — FIR低通滤波效果（时域）');
xlabel('时间/s'); ylabel('幅值'); grid on;
sgtitle('高斯白噪声滤除 — 低通FIR（Hamming窗）滤波器');

% Figure 6：FIR低通滤波前后频域对比
figure(6);
compare_spectra(Signal1, Sig1_FIR, N, f, Fs, ...
    '含高斯白噪声信号频谱', 'FIR低通滤波后频谱', ...
    '高斯白噪声滤除 — 低通FIR频域对比');

%% =========================================================================
%  第3节：工频干扰滤除 — 带阻滤波（IIR）
%  设计指标：fp = [45,55] Hz，fs = [48,52] Hz，Ap = 1 dB，As = 15 dB
%% =========================================================================
fp_bs = [45, 55];   % 带阻通带边界 (Hz)
fs_bs = [48, 52];   % 带阻阻带边界 (Hz)（包含50Hz工频）

Wp_bs = fp_bs / (Fs/2);    % 归一化通带截止频率
Ws_bs = fs_bs / (Fs/2);    % 归一化阻带截止频率

% 设计Butterworth带阻IIR滤波器
[ord_iir_bs, Wn_iir_bs] = buttord(Wp_bs, Ws_bs, Ap, As);
[b_iir_bs,   a_iir_bs ] = butter(ord_iir_bs, Wn_iir_bs, 'stop');

Sig2_IIR = filter(b_iir_bs, a_iir_bs, Signal2);

% Figure 7：IIR带阻滤波器幅频特性 + 滤波前后时域对比
figure(7);
subplot(2,1,1);
[h7, w7] = freqz(b_iir_bs, a_iir_bs, 2048, Fs);
plot(w7, 20*log10(abs(h7)+eps), 'k', 'LineWidth', 1.2);
title(sprintf('带阻IIR（Butterworth）幅频特性  阶数 N = %d', ord_iir_bs));
xlabel('频率/Hz'); ylabel('幅值/dB'); xlim([0 Fs/2]); grid on;
xline(fp_bs(1), '--b', sprintf('f_{p1}=%d Hz', fp_bs(1)));
xline(fp_bs(2), '--b', sprintf('f_{p2}=%d Hz', fp_bs(2)));
xline(fs_bs(1), '--r', sprintf('f_{s1}=%d Hz', fs_bs(1)));
xline(fs_bs(2), '--r', sprintf('f_{s2}=%d Hz', fs_bs(2)));

subplot(2,1,2);
plot(t0, Signal2, 'Color', [0.6 0.6 0.6]); hold on;
plot(t0, Sig2_IIR, 'k', 'LineWidth', 1.2); hold off;
legend('含工频干扰', 'IIR带阻滤波后', 'Location', 'best');
title('工频干扰 — IIR带阻滤波效果（时域）');
xlabel('时间/s'); ylabel('幅值'); grid on;
sgtitle('工频干扰滤除 — 带阻IIR（Butterworth）滤波器');

% Figure 8：IIR带阻滤波前后频域对比
figure(8);
compare_spectra(Signal2, Sig2_IIR, N, f, Fs, ...
    '含工频干扰信号频谱', 'IIR带阻滤波后频谱', ...
    '工频干扰滤除 — 带阻IIR频域对比');

%% =========================================================================
%  第4节：基线漂移滤除 — 高通滤波
%  设计指标：fp = 0.7 Hz，fs = 0.3 Hz，Ap = 1 dB，As = 15 dB
%
%  注：高通滤波器的通带截止频率 fp 大于阻带截止频率 fs，
%      即 fp=0.7Hz（高于此频率通过），fs=0.3Hz（低于此频率截止）
%% =========================================================================
fp_hp = 0.7;    % 高通通带截止频率 (Hz)
fs_hp = 0.3;    % 高通阻带截止频率 (Hz)

Wp_hp = fp_hp / (Fs/2);    % 归一化通带截止频率
Ws_hp = fs_hp / (Fs/2);    % 归一化阻带截止频率

%% ---- 4.1 高通IIR滤波器（Butterworth） ----
[ord_iir_hp, Wn_iir_hp] = buttord(Wp_hp, Ws_hp, Ap, As);
[b_iir_hp,   a_iir_hp ] = butter(ord_iir_hp, Wn_iir_hp, 'high');

Sig3_IIR = filter(b_iir_hp, a_iir_hp, Signal3);

% Figure 9：IIR高通滤波器幅频特性 + 滤波前后时域对比
figure(9);
subplot(2,1,1);
[h9, w9] = freqz(b_iir_hp, a_iir_hp, 2048, Fs);
plot(w9, 20*log10(abs(h9)+eps), 'k', 'LineWidth', 1.2);
title(sprintf('高通IIR（Butterworth）幅频特性  阶数 N = %d', ord_iir_hp));
xlabel('频率/Hz'); ylabel('幅值/dB'); xlim([0 Fs/2]); grid on;
xline(fp_hp, '--b', sprintf('f_p=%.1f Hz', fp_hp));
xline(fs_hp, '--r', sprintf('f_s=%.1f Hz', fs_hp));

subplot(2,1,2);
plot(t0, Signal3, 'Color', [0.6 0.6 0.6]); hold on;
plot(t0, Sig3_IIR, 'k', 'LineWidth', 1.2); hold off;
legend('含基线漂移', 'IIR高通滤波后', 'Location', 'best');
title('基线漂移 — IIR高通滤波效果（时域）');
xlabel('时间/s'); ylabel('幅值'); grid on;
sgtitle('基线漂移滤除 — 高通IIR（Butterworth）滤波器');

% Figure 10：IIR高通滤波前后频域对比
figure(10);
compare_spectra(Signal3, Sig3_IIR, N, f, Fs, ...
    '含基线漂移信号频谱', 'IIR高通滤波后频谱', ...
    '基线漂移滤除 — 高通IIR频域对比');

%% ---- 4.2 高通FIR滤波器（Hamming窗） ----
% 注：由于基线漂移过渡带极窄（Δf = 0.4 Hz），FIR高通滤波器阶数会非常高，
%     这是FIR滤波器的固有局限。实际工程中通常优先选用IIR高通滤波器。
delta_f_hp = (fp_hp - fs_hp) / Fs;     % 归一化过渡带宽（高通：fp > fs）
ord_fir_hp = ceil(6.6 / delta_f_hp);   % 估算FIR阶数（Hamming窗经验公式）

% 高通FIR滤波器要求阶数为偶数（保证对称型滤波器在奈奎斯特频率处有响应）
if mod(ord_fir_hp, 2) ~= 0
    ord_fir_hp = ord_fir_hp + 1;
end

fprintf('高通FIR滤波器阶数：%d（因过渡带极窄，阶数较高，运算时间稍长）\n', ord_fir_hp);

Wc_fir_hp = (fp_hp + fs_hp) / 2 / (Fs/2);   % 截止频率取中点
b_fir_hp  = fir1(ord_fir_hp, Wc_fir_hp, 'high', hamming(ord_fir_hp+1));

Sig3_FIR = fftfilt(b_fir_hp, Signal3);

% Figure 11：FIR高通滤波器幅频特性 + 滤波前后时域对比
figure(11);
subplot(2,1,1);
[h11, w11] = freqz(b_fir_hp, 1, 2048, Fs);
plot(w11, 20*log10(abs(h11)+eps), 'k', 'LineWidth', 1.2);
title(sprintf('高通FIR（Hamming窗）幅频特性  阶数 N = %d', ord_fir_hp));
xlabel('频率/Hz'); ylabel('幅值/dB'); xlim([0 Fs/2]); grid on;
xline(fp_hp, '--b', sprintf('f_p=%.1f Hz', fp_hp));
xline(fs_hp, '--r', sprintf('f_s=%.1f Hz', fs_hp));

subplot(2,1,2);
plot(t0, Signal3, 'Color', [0.6 0.6 0.6]); hold on;
plot(t0, Sig3_FIR, 'k', 'LineWidth', 1.2); hold off;
legend('含基线漂移', 'FIR高通滤波后', 'Location', 'best');
title('基线漂移 — FIR高通滤波效果（时域）');
xlabel('时间/s'); ylabel('幅值'); grid on;
sgtitle('基线漂移滤除 — 高通FIR（Hamming窗）滤波器');

% Figure 12：FIR高通滤波前后频域对比
figure(12);
compare_spectra(Signal3, Sig3_FIR, N, f, Fs, ...
    '含基线漂移信号频谱', 'FIR高通滤波后频谱', ...
    '基线漂移滤除 — 高通FIR频域对比');

fprintf('\nECG信号降噪处理完成，共生成12幅图像。\n');
fprintf('图1-2：信号生成与噪声对比\n');
fprintf('图3-6：高斯白噪声滤除（低通IIR + 低通FIR）\n');
fprintf('图7-8：工频干扰滤除（带阻IIR）\n');
fprintf('图9-12：基线漂移滤除（高通IIR + 高通FIR）\n');

%% =========================================================================
%%  辅助函数
%% =========================================================================

function plot_four_spectra(s0, s1, s2, s3, N, f)
% 绘制四路信号的单边幅度谱（频率范围 0~150 Hz）
    sigs   = {s0, s1, s2, s3};
    titles = {'干净ECG信号频谱', '加高斯白噪声信号频谱', ...
              '加工频干扰信号频谱', '加基线漂移信号频谱'};
    for k = 1:4
        mag = abs(fft(sigs{k})) / N;
        subplot(4,1,k);
        plot(f, mag, 'k');
        title(titles{k}); xlabel('频率/Hz'); ylabel('幅值');
        xlim([0 150]); grid on;
    end
end

function compare_spectra(sig_n, sig_f, N, f, Fs_val, title1, title2, fig_title)
% 绘制含噪信号与滤波后信号的频域幅度谱对比
    mag_n = abs(fft(sig_n)) / N;
    mag_f = abs(fft(sig_f)) / N;

    subplot(2,1,1);
    plot(f, mag_n, 'k');
    title(title1); xlabel('频率/Hz'); ylabel('幅值');
    xlim([0 Fs_val/2]); grid on;

    subplot(2,1,2);
    plot(f, mag_f, 'k');
    title(title2); xlabel('频率/Hz'); ylabel('幅值');
    xlim([0 Fs_val/2]); grid on;

    sgtitle(fig_title);
end

function signal = gen_ecg(N, Fs)
% 生成合成ECG信号（模拟真实心电图波形）
%   N  — 总采样点数
%   Fs — 采样率 (Hz)
%
% ECG波形由P波、QRS复合波（Q、R、S三个分量）和T波组成，
% 各分量均用高斯脉冲建模。
    HEART_RATE = 75;                           % 心率 (次/分钟)，正常静息范围 60-100
    signal     = zeros(1, N);
    period     = round(Fs * 60 / HEART_RATE);  % 每个心跳周期的采样点数

    for i = 0:floor(N/period)-1
        idx_start = i * period + 1;
        idx_end   = min(idx_start + period - 1, N);
        L         = idx_end - idx_start + 1;
        t_beat    = (0:L-1) / Fs;              % 当前心跳内的本地时间轴 (s)
        T         = period / Fs;               % 心跳周期 (s)

        % 各波形分量（高斯脉冲建模）
        beat = ...
            0.12 * gauss_pulse(t_beat, 0.10*T, 0.025*T) + ...  % P波
            1.00 * gauss_pulse(t_beat, 0.30*T, 0.008*T) - ...  % R波（QRS主峰）
            0.25 * gauss_pulse(t_beat, 0.27*T, 0.008*T) - ...  % Q波
            0.20 * gauss_pulse(t_beat, 0.33*T, 0.008*T) + ...  % S波
            0.35 * gauss_pulse(t_beat, 0.47*T, 0.040*T);       % T波

        signal(idx_start:idx_end) = signal(idx_start:idx_end) + beat;
    end
end

function y = gauss_pulse(t, t0, sigma)
% 生成以 t0 为中心、标准差为 sigma 的高斯脉冲
    y = exp(-((t - t0).^2) / (2 * sigma^2));
end

clear;
clc;

% 混凝土声速
v = 4000;

% 设定发射点为第10个阵元的坐标
A_x = 10 .* 2e-3;

% 可能的发射点和接收点的横坐标
% 阵元间隔为2mm，阵元长度为40mm
S_x = (0:20) .* 2e-3;

% 假设缺陷方向与水平方向平行
% 假设(R_start_x, R_i_y)为缺陷左端点，(R_end_x, R_i_y)为缺陷右端点
% 左端点的横坐标为0.015，右端点的横坐标为0.025
R_start_x = 15e-3;
R_end_x = 25e-3;
R_i_y = 50e-3;

% 发射波的参数
% 波频，单位为Hz
f0 = 1e6;

% 波形持续时间，单位为s
T = 5e-6;

% 步长
t_step = 1e-8;

% 时间序列
t = 0:t_step:200e-6;

% 发射波形
s = (1/2) * (heaviside(t) - heaviside(t - T)) .* (1 + cos(2 * pi / T * (t - T / 2))) .* cos(2 * pi * f0 * (t - T / 2));

% 用于储存水平缺陷的反射波形
wave = zeros(length(S_x), length(t));

% 用于储存反射波形的包络波形
wave_envelope = zeros(length(S_x), length(t));

% 储存声波经过路径的时间，也就是延时
t_delay = zeros(length(S_x), 1);

% 储存圆心横纵坐标C_x，C_y和半径d_i
circle = zeros(length(S_x), 3);

for i = 1:length(S_x)
    
	% 当前接收点的位置
    S_i_x = S_x(i);
    
	% 反射点坐标
    R_i_x = (A_x + S_i_x) / 2; 
	
    % 初始化d_i
    d_i = 0;
    
	% 裂纹面上的镜面反射
    if R_i_x >= R_start_x && R_i_x <= R_end_x
		% 计算声波经过的路径长度
		d_i = sqrt((R_i_x - A_x).^2 + R_i_y.^2) + sqrt((R_i_x - S_i_x).^2 + R_i_y.^2);
		% 将圆心坐标和 d_i 存储到数组中
        circle(i, 1) = S_i_x; 	% 圆心 x 坐标
        circle(i, 2) = 0; 		% 圆心 y 坐标
        circle(i, 3) = d_i;   	% 圆的半径d_i
	end
	
	% 计算声波经过路径的时间，也就是延时
	t_delay(i) = d_i / v;
	
	% 转换为时间步长的数量
    tn = round(t_delay(i) / t_step);

    % 增加延时，并存储到数组中
    wave(i, :) = wave(i, :) + [zeros(1, tn), s(1:length(s) - tn)];
	
    % 计算接收到的信号在时域中的包络
    wave_envelope(i, :) = wave_envelope(i, :) + abs(hilbert([zeros(1, tn), s(1:length(s) - tn)]));
end

%图3.2 发射波形
figure;
subplot(1,1,1);
plot(t, s);
xlabel('Time (s)');
ylabel('Amplitude');
title('Transmitted Waveform');
grid on;
xlim([0 15e-6]);
ylim([-1.1 1.1]);

%图3.3 换能器A自发自收
figure;
subplot(1,1,1);
hold on;
plot(t, s);
for i = 1:length(S_x)
    if i == 11
        plot(t, wave(i, :)); % 绘制每个接收点接收到的波形   
    end
end
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
title('Transmitted Waveform and Received Waveform at A');
grid on;
xlim([0 50e-6]);
ylim([-1.1 1.1]);
box on

%图3.4 全部换能器接收(水平)
figure;
subplot(1,1,1);
hold on;
for i = 1:length(S_x)
    if i >= 6 && i <=16
        plot (t, wave(i, :) ./ (2 * max(abs(wave(i, :))))+ i);
        %plot (t, wave_envelope(i, :) ./ (2 * max(abs(wave_envelope(i, :)))) + i);
        hold on;
    end
end
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
title('Received Waveform');
grid on;
xlim([25e-6 40e-6]);
ylim([5 17]);
box on
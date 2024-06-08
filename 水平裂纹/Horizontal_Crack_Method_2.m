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

% 存储圆的交点的数组
circle_intersection = [];

% 计算每两个圆的交点
for i = 1:length(S_x)
    for j = (i+1):length(S_x)
        % 获取第一个圆的信息
        C_1_x = circle(i, 1);
        C_1_y = circle(i, 2);
        C_1_r = circle(i, 3);
        
        % 获取第二个圆的信息
        C_2_x = circle(j, 1);
        C_2_y = circle(j, 2);
        C_2_r = circle(j, 3);
        
        % 计算圆心之间的距离
        l = sqrt((C_1_x - C_2_x).^2 + (C_1_y - C_2_y).^2);

        % 检查两个圆是否相交或相切
        if l > C_1_r + C_2_r || l < abs(C_1_r - C_2_r)
            % 没有交点
            continue;
        else
            % 计算圆的交点横坐标和纵坐标
            m = (C_1_r^2 - C_2_r^2 + l^2) / (2*l);
            n = sqrt(C_1_r^2 - m^2);
			
            x1 = C_1_x + m * (C_2_x - C_1_x) / l + n * (C_2_y - C_1_y) / l;
            y1 = C_1_y + m * (C_2_y - C_1_y) / l - n * (C_2_x - C_1_x) / l;
            x2 = C_1_x + m * (C_2_x - C_1_x) / l - n * (C_2_y - C_1_y) / l;
            y2 = C_1_y + m * (C_2_y - C_1_y) / l + n * (C_2_x - C_1_x) / l;

            % 判断交点是否在混凝土路面内部一侧
            if y1 > 0
                circle_intersection = [circle_intersection; [x1, y1]];
			end
			
            if y2 > 0
                circle_intersection = [circle_intersection; [x2, y2]];
			end
        end
    end
end
	
% 计算交点数量
num_circle_intersection = size(circle_intersection, 1);

% 初始化镜像点M坐标
M_x = 0;
M_y = 0;

% 累加交点的坐标
for i = 1:num_circle_intersection
    M_x = M_x + circle_intersection(i, 1);
    M_y = M_y + circle_intersection(i, 2);
end

% 计算交点的重心的坐标，也就是镜像点M坐标
M_x = M_x / num_circle_intersection;
M_y = M_y / num_circle_intersection;

% 计算D的坐标
D_x = (A_x + M_x)/2;
D_y = (0 + M_y)/2;

% 定义符号变量
syms x y; 

% 计算缺陷所在直线的方程直线一般式方程的参数，A、B、C 是直线的系数
%if S_i_x ~= M_x
k = (M_y - 0) / (M_x - A_x);
m = -(1 / k);
%end

%{
if S_i_x == M_x
    k = 0;
    m = D_y;
end
%}

% 直线方程
%line_eq = A * x + B * y + C;

% 储存缺陷所在直线方程和镜面反射椭圆的交点坐标
ellipse_intersection = [];

% 初始化 c 数组
a = zeros(length(S_x), 1);
c = zeros(length(S_x), 1);
b = zeros(length(S_x), 1);

% 椭圆方程
for i = 1:length(S_x)
    
	% 当前接收点的位置
    S_i_x = S_x(i);
    
	% 反射点坐标
    R_i_x = (A_x + S_i_x) / 2;	
        
        if R_i_x >= R_start_x && R_i_x <= R_end_x
            % 椭圆方程的参数
		
            % 椭圆中心点的坐标
            x_mid = R_i_x;
            y_mid = 0;

        
            % 椭圆的半长轴a，焦距c和半短轴b
            a(i) = t_delay(i) * v / 2;
            c(i) = abs(A_x - S_i_x) / 2;
            b(i) = sqrt(a(i).^2 - c(i).^2);

       
            % 椭圆方程
            %ellipse_eq = ((x - x_mid).^2 ./ a.^2) + ((y - y_mid).^2 ./ b.^2) - 1;
      
            % 求解直线和椭圆的交点
            % ellipse_intersection_point = solve([y - D_y == m * (x - D_y), ((x - x_mid).^2 ./ a.^2) + ((y - y_mid).^2 ./ b.^2) == 1], [x, y]);
            % 提取交点的坐标
            % x_intersection = double(ellipse_intersection_point.x);
            % y_intersection = double(ellipse_intersection_point.y);
            x_1 = ( x_mid ./ a(i) .^ 2 - m * (D_y - m * D_x) ./ b(i) .^ 2 + sqrt (1 ./ a(i) .^ 2 + m .^ 4 ./ b(i) .^ 2 - 1 ./ (a(i) .^ 2 .* b(i) .^ 2) .* (2 * m * x_mid + (D_y - m * D_x) .^ 2))) ./ (1 ./ a(i) .^2 + m .^ 2 ./ b(i) .^ 2);
            y_1 = m * (x_1 - D_x) + D_y;
            
            x_2 = ( x_mid ./ a(i) .^ 2 - m * (D_y - m * D_x) ./ b(i) .^ 2 - sqrt (1 ./ a(i) .^ 2 + m .^ 4 ./ b(i) .^ 2 - 1 ./ (a(i) .^ 2 .* b(i) .^ 2) .* (2 * m * x_mid + (D_y - m * D_x) .^ 2))) ./ (1 ./ a(i) .^2 + m .^ 2 ./ b(i) .^ 2);
            y_2 = m * (x_2 - D_x) + D_y;
            
        
            % 将交点坐标存储到二维数组中
            ellipse_intersection = [ellipse_intersection; x_1, y_1];
            ellipse_intersection = [ellipse_intersection; x_2, y_2];     
	end
end

figure
hold on;
scatter(ellipse_intersection(:,1), ellipse_intersection(:,2), 150, 'k', 'x');
plot(ellipse_intersection(:,1), ellipse_intersection(:,2),'b');
hold off;
xlabel('X (m)');
ylabel('Y (m)');
title('Intersection and Crack Line');
legend('Intersection', 'Crack Line');
xlim([0.014 0.026]);
ylim([0.045 0.055]);
set(gca, 'YDir', 'reverse');
box on;

% 绘制直线方程和所有椭圆方程的图像
figure;
hold on;

% 绘制直线方程
ezplot(y - D_y == -(1/k) * (x - D_y), [-0.1 0.1 -0.1 0.1]); % 根据需要调整坐标轴范围

% 绘制椭圆方程
for i = 1:length(S_x)
    % 当前接收点的位置
    S_i_x = S_x(i);
    
    % 反射点坐标
    R_i_x = (A_x + S_i_x) / 2;
	
    if R_i_x >= R_start_x && R_i_x <= R_end_x
        % 椭圆方程的参数
        
        % 椭圆中心点的坐标
        x_mid = R_i_x;
        y_mid = 0;

        % 椭圆的半长轴a，焦距c和半短轴b
        a(i) = t_delay(i) * v / 2;
        c(i) = abs(A_x - S_i_x) / 2;
        b(i) = sqrt(a(i).^2 - c(i).^2);

        % 绘制当前椭圆的图形
        ezplot(@(x,y) ((x - x_mid).^2 / a(i).^2) + ((y - y_mid).^2 / b(i).^2) - 1, [-0.1 0.1 -0.1 0.1]); % 根据需要调整坐标轴范围
        
		% 反向Y坐标轴的方向
        set(gca, 'YDir', 'reverse');
    end
end

title('Intersection of Line and Ellipse');
xlabel('X');
ylabel('Y');
legend('Line', 'Ellipse'); % 添加图例
box on;
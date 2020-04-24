%% 微信公众号：数学模型（MATHmodels）
%  联系方式：mathmodels@outlook.com

%% 不倒翁半球数据
r = 0.5;             % [m] 半径
m = 1.0              % [kg] 质量
g = 9.8;             % [m/s] 重力加速度
h = 3/8*r;           % [m] 重心位置
Ic = 83/320*m*r^2;   % [kg*m^2] 转动惯量

%% 仿真参数
tmax = 2.3*5;           % [s] 仿真时间
theta0 = deg2rad(60);% [deg] 初始角度
omega0 = 0;          % [deg/s] 初始角速度
xo = -theta0*r;      % [m] 初始位移

%% 画不倒翁，为了更像玩具不倒底，加了个三角形的顶
Xb = r*[cosd(180:360) 0];
Yb = r*[sind(180:360) 2];
[xb, yb] = rotxy(Xb,Yb, theta0);
hb = fill(xb+xo,yb,[1,0.8,0.8]);
axis image; axis([-2.25,2.25,-0.5,2]);

%% 调用 ode45 求解微分方程
[t, to] =ode45(@odes,[0,tmax], [theta0 omega0], [], r, m, g, h, Ic);
title(sprintf('t = %8.4f', 0));

for i =2:length(t)
    theta = to(i,1);
    
    % 更新中心位移
    xo = -theta*r;
    
    % 根据角度旋转座标点
    [xb, yb] = rotxy(Xb, Yb, theta);
    
    % 更新图形并显示
    set(hb, 'XData',xb+xo, 'YData',yb)
    title(sprintf('t = %8.4f', t(i)));
    pause(t(i)-t(i-1)) % 暂停 dt 秒
end

% -------------------------------------------------------------------------

function dy = odes(t, y, r, m, g, h, Ic)
% y(1) = theta; y(2) = d(theta)/dt

nume = -h*m*(g+r*y(2)^2)*sin(y(1));    % 分子
deno = Ic+m*(r^2+h^2-2*h*r*cos(y(1))); % 分母

dy = [y(2); nume/deno];
end

% -------------------------------------------------------------------------

function [x, y] = rotxy(x0, y0, theta)
% 将 (x0,y0) 绕 (0,0) 旋转 theta 度
x = x0*cos(theta) - y0*sin(theta);
y = x0*sin(theta) + y0*cos(theta);
end
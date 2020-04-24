%% 微信公众号：数学模型（MATHmodels）
%  联系方式：mathmodels@outlook.com

clear; clc
%% 可调参数
T = 2.18;           % [s]      小姐姐摇晃周期
amp = deg2rad(30);  % [rad]    小姐姐摇晃幅度
a0 = 0.01;          % [m]      滚动阻力系数
dt = 1e-2;          % [s]      仿真步长

%% 物理常数
g = 9.8;            % [m/s]    重力加速度

%% 底座（含小姐姐的下半身）数据
r = 0.5;            % [m]      底座半径
m1 = 275;           % [kg]     底座（含小姐姐的下半身）质量
h = 0.2442;         % [m]      底座（含小姐姐的下半身）重心位置
I1 = 28.34;         % [kg*m^2] 底座（含小姐姐的下半身）转动惯量
theta = deg2rad(0); % [rad]    初始角度
w1 = 0;             % [deg/s]  初始角速度
H = 0.9;            % [m]      小姐姐下半身长度

%% 小姐姐上半身数据
m2 = 20;            % [kg]     小姐姐上半身质量
l = 0.73/2;         % [m]      小姐姐上半身长度的一半（质心到腰的距离）
I2 = 0.94;          % [kg*m^2] 小姐姐上半身转动惯量

[beta, w2, dw2] = sway(0, T, amp);
hp = plotumbler(theta, beta);

t = 0:dt:12*3600;

for i = 2:length(t)
    [beta(i), w2, dw2] = sway(t(i), T, amp);
    
    theta(i) = theta(i-1) + w1 * dt;
    
    dw1 = odes(theta(i),w1, beta(i),w2,dw2,a0,m1,m2,H,h,r,l,I1,I2,g);
    
    w1 = w1 + dw1 * dt;
    
    [F(i,:), M(i,:)] = force(m2,r,H,l,I2,g,w1,w2,dw1,dw2,theta(i),beta(i));
    
    plotumbler(theta(i), beta(i), hp);
end

% -------------------------------------------------------------------------

function [F, M] = force(m2, r, H, l, I2, g, w1, w2, dw1, dw2, theta, beta)
% 求解小姐姐腰部作用力和力偶
Fx = -m2*(r*dw1 + l*cos(beta + theta)*(dw2 + dw1) ...
          - H*sin(theta)*w1^2 + H*cos(theta)*dw1 ...
          - l*sin(beta + theta)*(w2 + w1)^2        );
      
Fy = -m2*(   -g + l*sin(beta + theta)*(dw2 + dw1) ...
          + H*cos(theta)*w1^2 + H*sin(theta)*dw1 ...
          + l*cos(beta + theta)*(w2 + w1)^2        );

F = [Fx, Fy];
M = l*cos(theta+beta)*Fx + l*sin(theta+beta)*Fy - I2*(dw1+dw2);
end

% -------------------------------------------------------------------------

function dw1 = odes(theta,w1, beta,w2,dw2, a0, m1,m2, H,h,r,l, I1,I2, g)
a = a0 * sign(w1); % d = d0 * w1;
num = l*m2*(g+r*(w1+w2)^2+dw2*a)*sin(beta + theta) ...
      + l*m2*(  a*(w1+w2)^2-dw2*r)*cos(beta + theta) ...
      - (m1+m2)*a*g - dw2*l^2*m2 - I2*dw2 - (m1+m2)*g*h*sin(theta) ...
      + dw2*h*l*m2*cos(beta) - l*m2*(h*(w1+w2)^2+H*w1^2)*sin(beta) ...
      + w1^2*(H*m2 - m1*h)*(a*cos(theta)+ r*sin(theta));

den = I1 + I2 + m1*(h^2 + r^2) + m2*(l^2 + r^2 - H*h) ...
      + 2*l*m2*r*cos(beta + theta) - a*l*m2*sin(beta + theta) ...
      + (H-h)*l*m2*cos(beta) + (H*m2 -h*m2- 2*m1*h)*r*cos(theta) ...
      + (m1*h-H*m2)*a*sin(theta);
      
dw1 = num/den;
end
% -------------------------------------------------------------------------

function [beta, w, dw] = sway(t, T, amp)
if nargin==1; T = 2; amp = deg2rad(30); end
omega = 2*pi/T;                % [rad/s] 摇晃周期的角速度
beta = amp*sin(omega*t);       % [rad]   摇晃角度
w = amp*omega*cos(omega*t);    % = d(beta)/dt 摇晃角速度 
dw = -amp*omega^2*sin(omega*t);% = d(w)/dt 摇晃角加速度
end

% -------------------------------------------------------------------------

function h = plotumbler(theta, beta, h)
if nargin==0; theta = deg2rad(60); beta = deg2rad(-30); end
r = 0.5;                       % [m] 半径
xo = -theta* r;                % [m] 底座半球中心位置
%　底座半球坐标
xb = r*cosd(180:360); yb = r*sind(180:360);
[xb, yb] = rotxy(xb,yb, theta);
xb = xb + xo;
% 小姐姐坐标
x = [0.0  0.6  0.3  0.2  0.0 -0.2 -0.3 -0.6 ...
     0.1  0.1  0.1  0.0 -0.1 -0.1 -0.1  NaN];
y = [1.6  1.2  1.2  1.4  1.4  1.4  1.2  1.2 ...
     0.0  0.5  0.9  1.0  0.9  0.5  0.0  NaN];
 
I = [1 5 12 16 2 3 4 5 6 7 8 16 9 10 11 12 13 14 15];
It = 1:8; % 表示上半身节点编号
Ic = 12;  % 表示腰的节点（小姐姐摇晃的转动中心）编号

[x, y] = rotxy(x,y,theta); % 小姐姐整体随底座转动 theta 角
% 小姐姐上半身绕腰转动 beta 角
[x(It), y(It)] = rotxy(x(It), y(It), beta, x(Ic), y(Ic));
xm = x(I)+xo; ym = y(I);
if nargin<3
    h(1) = fill(xb, yb, [1,0.8,0.8]); hold on;
    h(2) = plot(xm,ym,'bo-','linewidth',2, 'markersize',10);
    axis image; axis([-2.25,2.25,-0.5,2]);
else
    set(h(1), 'XData',xb, 'YData',yb)
    set(h(2), 'XData',xm, 'YData',ym)
end
drawnow
end

% -------------------------------------------------------------------------

function [x, y] = rotxy(x0, y0, theta, xc, yc)
% 将 (x0,y0) 绕 (xc,yc) 旋转 theta 
if nargin==3; xc = 0; yc = 0; end
x0 = x0 - xc; y0 = y0 - yc;
x =  x0*cos(theta) - y0*sin(theta) + xc;
y =  x0*sin(theta) + y0*cos(theta) + yc;
end

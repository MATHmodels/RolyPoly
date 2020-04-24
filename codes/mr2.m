%% 微信公众号：数学模型（MATHmodels）
%  联系方式：mathmodels@outlook.com

[hc, Ic, mc] = sphcapci(); % 底座球壳质心、转动量和质量
[hw, Iw, mw] = sphshlci(); % 底座配重质心、转动量和质量
[hh, Ih, mh] = cylci();    % 小姐姐下半身质心、转动量和质量

m = mc + mw + mh;
h = (hc*mc + hw*mw + hh*mh)/m;
I = Ic + mc*(hc-h)^2 + Iw + mw*(hw-h)^2 + Ih + mh*(hh-h)^2;

% -------------------------------------------------------------------------

function [Zc, Ic, m] = sphshlci(R, r, m)
% 半球壳的质心、转动量和质量

if nargin==0
    R = 0.5;
    r = R - 0.005;
    rho = 7.87e3;  % kg/m^3
    m = rho*2/3*pi*(R.^3 - r.^3);
end

dr3 = (R.^3 - r.^3);
dr4 = (R.^4 - r.^4);
dr5 = (R.^5 - r.^5);

V = 2/3*pi*dr3;
Zc = -3/8 * dr4./dr3;
Iy = 2/5 * m * dr5/dr3;
Ic = Iy - m*Zc^2;
end

% -------------------------------------------------------------------------

function [Zc, Ic, m] = sphcapci(r, b, m)
% 球缺的质心、转动动量和质量

% syms r z b m
% dv = pi*(r^2-z^2); 
% V = int(dv, z, -r, -b)
% Zc = simplify(int(dv*z, z, -r, -b)/V)
% dm = m/V * dv;
% dI = 1/4*(r^2-z^2)*dm + z^2*dm;
% I = int(dI, z,-r,-b)

if nargin==0
    r = 0.495; b = -r/2; m = 188.81; 
end

Zc = 3*(b - r)^2/(4*b - 8*r);
Iy = (9*b^2 + 17*r^2 + 18*r^3/(b-2*r))*m/20;
Ic = Iy - m*Zc^2;
V = pi*(b + r)^2*(2*r-b)/3;
end

% -------------------------------------------------------------------------

function [Zc, Ic, m] = cylci(r, h, m)
% 圆柱的质心、转动动量和质量

if nargin==0;
    h = 0.9; % h = 0.73;
    r = 0.094; m = 45*h/1.63;
end

Zc = h/2;
Iy = 1/4*m*r^2 + 1/3*m*h^2;
Ic = Iy - m*Zc^2;
end

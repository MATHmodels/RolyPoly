%% 微信公众号：数学模型（MATHmodels）
%  联系方式：mathmodels@outlook.com

syms theta(t) r h F N g m I

x = -r*theta;
xo = [x; 0];
xc = xo + h*[sin(theta); -cos(theta)];
ac = diff(xc, t, 2);

eqF = [F; N-m*g]==m*ac;
eqI = -(h*sin(theta))*N + (r-h*cos(theta))*F == I*diff(theta,t,2);

[Fv,Nv] = solve(eqF, F, N)

eqI = subs(eqI, F, Fv);
eqI = subs(eqI, N, Nv);

ddt = isolate(eqI, diff(theta,t,2))

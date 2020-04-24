%% 微信公众号：数学模型（MATHmodels）
%  联系方式：mathmodels@outlook.com

syms theta(t) beta(t) r h F N g m1 m2 I1 I2 d H Fx Fy T l w1 w2 dw2

x = -r*theta;

xo = [x; 0];
xc = xo + h*[ sin(theta); -cos(theta)];
xd = xo + H*[-sin(theta);  cos(theta)];
xe = xd + l*[-sin(theta+beta); cos(theta+beta)];

ac = diff(xc, t, 2);
ad = diff(xd, t, 2);
ae = diff(xe, t, 2);


eqF1 = [F-Fx; N-m1*g-Fy]==m1*ac;
eqF2 = [  Fx; Fy-m2*g  ]==m2*ae;

eqI1 = -(h*sin(theta)+d)*N + (r-h*cos(theta))*F + M == I1*diff(theta,t,2);
eqI2 = l*cos(theta+beta)*Fx + l*sin(theta+beta)*Fy - M == I2*diff(theta+beta,t,2);

eqI = eqI1 + eqI2;

[Fv,Nv,Fxv, Fyv] = solve([eqF1; eqF2], F, N, Fx, Fy);

eqI = subs(eqI, {F, N, Fx, Fy}, {Fv, Nv, Fxv, Fyv});


ddt = isolate(eqI, diff(theta,t,2));
ddt = subs(ddt, diff(theta,t), w1);
ddt = subs(ddt, diff(beta,t), w2);
ddt = subs(ddt, diff(beta,t,t), dw2);

ddt = subs(ddt, {w2, dw2, d, I2, m2}, {0,0,0,0,0});

[num, den] = numden(rhs(ddt))

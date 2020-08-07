dimen=input('x的维度：');
tic
x0 = linspace(0,0,dimen);
% 负梯度法
[xk,yk,k] = minus_gradient(x0,dimen);
%非线性共轭梯度法
% [xk,yk,k] = con_gradient(x0,dimen);
% DFP方法
% [xk,yk,k] = DFP(x0,dimen);
runtime = toc;
dimen=input('x��ά�ȣ�');
tic
x0 = linspace(0,0,dimen);
% ���ݶȷ�
[xk,yk,k] = minus_gradient(x0,dimen);
%�����Թ����ݶȷ�
% [xk,yk,k] = con_gradient(x0,dimen);
% DFP����
% [xk,yk,k] = DFP(x0,dimen);
runtime = toc;
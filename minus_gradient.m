function [x_value,y_value,k] = minus_gradient(x0,n)
    %x0以及所有算法原本的列向量在这里均为行向量
    epsilon = 1e-5;
    xk = x0; k = 1;
    grad = -gradfunc(xk,n);
    while(norm(grad) >= epsilon)
        xk = Armijo(-grad,n,xk,grad);
%         alpha = wolfe(xk,n,grad);
%         xk = xk + alpha * grad;
        grad = -gradfunc(xk,n);
        k = k + 1;
    end
    x_value = xk;
    y_value = fun(xk,n);
end

% function [alpha] = wolfe(xk,n,dk)
%     rho = 0.25; sigma = 0.75;
% %     alpha = 1; a = 0; b = 1;
%     mk = 25; m = 1;
%     alpha = 0.8; delta = 0.8;
%     gxk = gradfunc(xk,n);
%     dkc = dk';
%     while (m<=mk)
%         if fun(xk+alpha*dk,n) <= fun(xk,n)+rho*alpha*gxk*dkc
%             if gradfunc(xk+alpha*dk,n)*dkc >= sigma*gxk*dkc
%                 break;
% %             else
% %                 a = alpha;
%             end        
% %         else
% %             b = alpha; 
%         end
% %         alpha =(a+b)/2;
%         alpha = alpha * delta;
%     end
% end
function newxk = Armijo(gfun_xk,n,xk,d)
%ARMIJO 此处显示有关此函数的摘要
%   此处显示详细说明
%gfun_xk,xk,d均为行向量
    mk = 1; max_mk = 20;
    rho = 0.5; sigma = 1e-6;
    while(mk <= max_mk)
        newxk = xk + rho^mk * d;
        if (fun(newxk,n) <= fun(xk,n) + sigma * rho^mk * gfun_xk * d')
            break;
        end
        mk = mk + 1;
    end
end

% function c = fun(x,n)
%     c = x(1)^2 - x(1) + 4*x(2)^2 +n-n;
% end
% 
% function c = gradfunc(x,n)
%     c = [0 0];
%     c(1) = 2*x(1) -1;
%     c(2) = 8*x(2);
% end

%返回数值
function yx = fun(x,n)
    yx = 0;
    for i=1:29
        ti = i/29;
        yx = yx + rix(x,ti,n)^2;
    end
    yx = yx + x(1)^2 + (x(2) - x(1)^2 - 1)^2;
end
%返回行向量
function grad = gradfunc(x,n)
    grad = linspace(0,0,n);
    rixc = linspace(0,0,n);
    mixc = linspace(0,0,n);
    tic = linspace(0,0,n);
    for i=1:29
        tic(i) = i/29;
        rixc(i) = rix(x,tic(i),n);
        mixc(i) = mix(x,tic(i),n);
    end
    for j=1:n
        for i=1:29
            grad(j) = grad(j) + 2*rixc(i)*((j-1)*tic(i)^(j-2)-2*mixc(i)*tic(i)^(j-1));
        end
    end
    grad(1) = grad(1) + 2*x(1) + 4*(x(1)^3 + x(1) - x(1)*x(2));
    grad(2) = grad(2) + 2 * (x(2) - x(1)^2 - 1);
end
%返回数值
function  yx = rix(x,ti,n)
    yx = 0;
    for j=2:n
        yx = yx + (j-1) * x(j) * ti^(j-2);
    end
    yx = yx - mix(x,ti,n)^2 - 1;
end
%返回数值
function yx = mix(x,ti,n)
    yx = 0;
    for j=1:n
        yx = yx + ti^(j-1) * x(j);
    end
end
function [x_value,y_value,k] = con_gradient(x0,n)
%CON_GRADIENT 此处显示有关此函数的摘要
%   此处显示详细说明
%x,x0均为行向量
%该函数内部均为行向量，可在返回后变为列向量
%     gfun = jacobian(fun,x);
    grad2 = -gradfunc(x0,n);
    epsilon = 1e-5;
    xk = x0; k = 1;
    %步长
    d = grad2;
    while(norm(grad2) >= epsilon)
        %xk梯度
        grad1 = grad2;
%         alpha = wolfe(xk,n,d);
%         xk = xk + alpha * d;
        xk = Armijo(-grad1,n,xk,d);
        %x(k+1)梯度
        grad2 = gradfunc(xk,n);
        d = -grad2 + norm(grad2)^2/norm(grad1)^2 * d;
        k = k+1;
    end
    x_value = xk;
    y_value = fun(xk,n);
end

% function [alpha] = wolfe(xk,n,dk)
%     rho = 0.25; sigma = 0.75;
%     alpha = 1; a = 0; b = Inf;
%     gxk = gradfunc(xk,n);
%     dkc = dk';
%     while (1)
%         if ~(fun(xk+alpha*dk,n)<=fun(xk,n)+rho*alpha*gxk*dkc)
%             b = alpha;
%             alpha = (alpha+a)/2;
%             continue;
%         end
%         if ~(gradfunc(xk+alpha*dk,n)*dkc >= sigma*gxk*dkc)
%             a = alpha;
%             alpha = min([2*alpha, (b+alpha)/2]);
%             continue;
%         end
%         break;
%     end
% end
function newxk = Armijo(gfun_xk,n,xk,d)
%ARMIJO 此处显示有关此函数的摘要
%   此处显示详细说明
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

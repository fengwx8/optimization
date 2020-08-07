function [x_value,y_value,k] = DFP(x0,n)
%DFP 此处显示有关此函数的摘要
%   此处显示详细说明
    epsilon = 1e-5;
    xk = x0; k = 1;
    g2 = gradfunc(xk,n);
    D = eye(n);
    while(norm(g2) >= epsilon)
        g1 = g2;
        d = -D * g1';
        lambda = Armijo(g1,n,xk,d');
        sk = lambda * d';
%         alpha = wolfe(xk,n,d');
%         sk = alpha * d';
        xk = xk + sk;
        g2 = gradfunc(xk,n);
        yk = g2 - g1;
        D = D + (sk'*sk)/(sk*yk') - (D*yk')*(yk*D)/(yk*D*yk');
        k = k + 1;
    end
    x_value = xk;
    y_value = fun(xk,n);
end


function [alpha] = wolfe(xk, n,dk)
    rho = 0.25; sigma = 0.75;
%     alpha = 1; a = 0; b = 1; 
    dkc = dk';
    gxk = gradfunc(xk,n);
    fxk = fun(xk,n);
    alpha = 0.5; delta = 0.5; 
    m = 1; mk = 20;
%     while b-a<1e-4
    while(m<=mk)
        if fun(xk + alpha*dk,n) <= fxk + rho*alpha*gxk*dkc
            if gradfunc(xk + alpha*dk,n)*dkc >= sigma*gxk*dkc
                break;
%             else
%                 a = alpha;
            end
%         else
%             b = alpha;
        end
%         alpha = (b+a)/2;
        m = m+1;
        alpha = alpha * delta;
    end

end
function lambda = Armijo(gfun_xk,n,xk,d)
%ARMIJO 此处显示有关此函数的摘要
%   此处显示详细说明
    mk = 0; max_mk = 20;
    rho = 0.5; sigma = 1e-6;
    while(mk <= max_mk)
        newxk = xk + rho^mk * d;
        if (fun(newxk,n) <= fun(xk,n) + sigma * rho^mk * gfun_xk * d')
            break;
        end
        mk = mk + 1;
    end
    lambda = rho^mk;
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

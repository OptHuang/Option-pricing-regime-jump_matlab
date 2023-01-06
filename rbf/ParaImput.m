function [T,K,sigma1,sigma2,r1,r2,lam1,lam2,L,epsilon,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,p,q,a1,a2,kappaM,kappaK] = ParaImput()
    % 由于参数比较多，所以写一个传入参数的函数方便修改调用
    
    T = 1;  % 传入时间上界T
    K = 10;  % 传入K值 
    
    % 传入两种制度下的波动率sigma,无风险利率r,红利d,跳跃强度λ,期望k
    sigma1 = 0.3;
    sigma2 = 0.4;
    r1 = 0.05;
    r2 = 0.05;
    lam1 = 0;
    lam2 = 0;

    
    % 截取边界
    L = log(3*K);

    epsilon = 1e-6;
    
    % 传入不同制度之间的转化矩阵A
    a11 = -3;
    a12 = 3;
    a21 = 2;
    a22 = -2;
    A = [a11 a12;a21 a22];
    
    x0 = -L;  % 传入空间左侧起点
    t0 = 0;   % 传入时间起点
    
    Nt = 200;  % 传入时间分划段数
    Nx = 300;   % 传入空间分划段数

    dt = T/Nt;    % 计算出时间分划间距
    dx = 2*L/Nx;  % 计算出空间分划间距
    
    % Merton型的参数
    mu = -0.025;
    delt = sqrt(0.05);  
    kappaM = exp(mu+delt^2/2)-1;
    
    % Kou型参数在h函数内部
    p = 0.3445;  
    q = 1-p;
    a1 = 3.0465;
    a2 = 3.0775;
    kappaK = p*a1/(a1-1)+(1-p)*a2/(a2+1)-1;
    
end
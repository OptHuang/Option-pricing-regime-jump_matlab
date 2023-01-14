function [T,K,sigma1,sigma2,r1,r2,lam1,lam2,gamma,mun,rho0,eps,epsilon,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2,p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = ParaImput()
    % 由于参数比较多，所以写一个传入参数的函数方便修改调用
    
    T = 1;  % 传入时间上界T
    K = 1;  % 传入K值 
    
    % 传入两种制度下的波动率sigma,无风险利率r,红利d,跳跃强度λ,期望k
    sigma1 = 0.8;
    sigma2 = 0.3;
    r1 = 0.05;
    r2 = 0.05;
    d1 = 0;
    d2 = 0;
    lam1 = 0.2;
    lam2 = 0.2;

    
    % 传入pcm的几个参数
    gamma = 0.9;
    mun = 0.4;
    rho0 = 1.5;
    eps = 1e-9;
    
    % RBF 参数
    epsilon = 1e-6;
    
    % 传入不同制度之间的转化矩阵A
    a11 = -2;
    a12 = 2;
    a21 = 3;
    a22 = -3;
    A = [a11 a12;a21 a22];
    
    % 截取边界
    L = TruncationTech(K,0,sigma1,sigma2,r1,r2,d1,d2,lam1,lam2,A);
    
    x0 = -L;  % 传入空间左侧起点
    t0 = 0;   % 传入时间起点
    
    Nt = 500;  % 传入时间分划段数
    Nx = 256;   % 传入空间分划段数

    dt = T/Nt;    % 计算出时间分划间距
    dx = 2*L/Nx;  % 计算出空间分划间距
    
    % Merton型的参数
    [mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2] = Mertonpara(sigma1,sigma2,r1,r2,d1,d2,lam1,lam2,A);
    
    % Kou型参数在h函数内部
    [p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = Koupara(sigma1,sigma2,r1,r2,d1,d2,lam1,lam2,A);
    
end
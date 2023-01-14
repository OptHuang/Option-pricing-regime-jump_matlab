function [T,K,sigma1,sigma2,r1,r2,lam1,lam2,gamma,mun,rho0,eps,epsilon,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2,p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = ParaImput()
    % ���ڲ����Ƚ϶࣬����дһ����������ĺ��������޸ĵ���
    
    T = 1;  % ����ʱ���Ͻ�T
    K = 1;  % ����Kֵ 
    
    % ���������ƶ��µĲ�����sigma,�޷�������r,����d,��Ծǿ�Ȧ�,����k
    sigma1 = 0.8;
    sigma2 = 0.3;
    r1 = 0.05;
    r2 = 0.05;
    d1 = 0;
    d2 = 0;
    lam1 = 0.2;
    lam2 = 0.2;

    
    % ����pcm�ļ�������
    gamma = 0.9;
    mun = 0.4;
    rho0 = 1.5;
    eps = 1e-9;
    
    % RBF ����
    epsilon = 1e-6;
    
    % ���벻ͬ�ƶ�֮���ת������A
    a11 = -2;
    a12 = 2;
    a21 = 3;
    a22 = -3;
    A = [a11 a12;a21 a22];
    
    % ��ȡ�߽�
    L = TruncationTech(K,0,sigma1,sigma2,r1,r2,d1,d2,lam1,lam2,A);
    
    x0 = -L;  % ����ռ�������
    t0 = 0;   % ����ʱ�����
    
    Nt = 500;  % ����ʱ��ֻ�����
    Nx = 256;   % ����ռ�ֻ�����

    dt = T/Nt;    % �����ʱ��ֻ����
    dx = 2*L/Nx;  % ������ռ�ֻ����
    
    % Merton�͵Ĳ���
    [mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2] = Mertonpara(sigma1,sigma2,r1,r2,d1,d2,lam1,lam2,A);
    
    % Kou�Ͳ�����h�����ڲ�
    [p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = Koupara(sigma1,sigma2,r1,r2,d1,d2,lam1,lam2,A);
    
end
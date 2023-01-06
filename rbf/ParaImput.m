function [T,K,sigma1,sigma2,r1,r2,lam1,lam2,L,epsilon,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,p,q,a1,a2,kappaM,kappaK] = ParaImput()
    % ���ڲ����Ƚ϶࣬����дһ����������ĺ��������޸ĵ���
    
    T = 1;  % ����ʱ���Ͻ�T
    K = 10;  % ����Kֵ 
    
    % ���������ƶ��µĲ�����sigma,�޷�������r,����d,��Ծǿ�Ȧ�,����k
    sigma1 = 0.3;
    sigma2 = 0.4;
    r1 = 0.05;
    r2 = 0.05;
    lam1 = 0;
    lam2 = 0;

    
    % ��ȡ�߽�
    L = log(3*K);

    epsilon = 1e-6;
    
    % ���벻ͬ�ƶ�֮���ת������A
    a11 = -3;
    a12 = 3;
    a21 = 2;
    a22 = -2;
    A = [a11 a12;a21 a22];
    
    x0 = -L;  % ����ռ�������
    t0 = 0;   % ����ʱ�����
    
    Nt = 200;  % ����ʱ��ֻ�����
    Nx = 300;   % ����ռ�ֻ�����

    dt = T/Nt;    % �����ʱ��ֻ����
    dx = 2*L/Nx;  % ������ռ�ֻ����
    
    % Merton�͵Ĳ���
    mu = -0.025;
    delt = sqrt(0.05);  
    kappaM = exp(mu+delt^2/2)-1;
    
    % Kou�Ͳ�����h�����ڲ�
    p = 0.3445;  
    q = 1-p;
    a1 = 3.0465;
    a2 = 3.0775;
    kappaK = p*a1/(a1-1)+(1-p)*a2/(a2+1)-1;
    
end
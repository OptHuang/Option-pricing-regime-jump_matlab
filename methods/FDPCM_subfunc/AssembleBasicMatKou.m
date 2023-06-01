function [B, d1t0, d2t0, F0t0, Phi_star_t0] = AssembleBasicMatKou(dt,dx,Nx,sigma1,sigma2,nK1,nK2,mK1,mK2,A,lam1,lam2,L,K,p,q,a1,a2)

    % beta是时间空间剖分网比
    beta = dt/(dx^2);
    
    %%%%%%%%%%%%%%%%%形成矩阵B
    % 形成矩阵B1,B2
    B1 = zeros(Nx-1);
    B2 = zeros(Nx-1);
    for i = 2:Nx-2
        B1(i,i-1) = -beta*sigma1^2/2;
        B1(i,i) = 1+beta*sigma1^2;
        B1(i,i+1) = -beta*sigma1^2/2;
        B2(i,i-1) = -beta*sigma2^2/2;
        B2(i,i) = 1+beta*sigma2^2;
        B2(i,i+1) = -beta*sigma2^2/2;
    end
    B1(1,1) = 1+beta*sigma1^2;
    B1(1,2) = -beta*sigma1^2/2;
    B1(Nx-1,Nx-2) = -beta*sigma1^2/2;
    B1(Nx-1,Nx-1) = 1+beta*sigma1^2;
    B2(1,1) = 1+beta*sigma2^2;
    B2(1,2) = -beta*sigma2^2/2;
    B2(Nx-1,Nx-2) = -beta*sigma2^2/2;
    B2(Nx-1,Nx-1) = 1+beta*sigma2^2;
    
    % 组装矩阵B
    B = [B1 zeros(size(B1)); zeros(size(B1)) B2];


    %%%%%%%%%%%%%%%%%形成对角Di矩阵的对角元素组成的向量dit0

    z1 = ones(Nx-1,1);
    z2 = ones(Nx-1,1);
    for i = 1:Nx-1
        z1(i) = exp(i*dx*(nK2-nK1));
        z2(i) = 1/z1(i);
    end
    d1t0 = -A(1,2)*dt*exp(-(nK2-nK1)*L)*z1;
    d2t0 = -A(2,1)*dt*exp(-(nK1-nK2)*L)*z2;


    %%%%%%%%%%%%%%%%形成F0t0

    % 形成向量F0t0, F0t0=FDt0 -dt*lam*FIt0 -dt*Fct0
    % 形成FDt0
    w_cur1 = (K-exp(-L))*exp(-mK1*1*dt+nK1*L);
    w_cur2 = (K-exp(-L))*exp(-mK2*1*dt+nK2*L);
    FDt0 = zeros(2*Nx-2,1);
    FDt0(1) = -beta*sigma1^2/2*w_cur1;     
    FDt0(Nx) = -beta*sigma2^2/2*w_cur2;
    % 形成FIt0
    w_pre1 = (K-exp(-L))*exp(nK1*L);
    w_pre2 = (K-exp(-L))*exp(nK2*L);
    FIt0 = zeros(2*Nx-2,1);
    for i = 1:Nx-1
        FIt0(i) = lam1*h(0,i,0,nK1,dx,p,q,a1,a2)*w_pre1;
        FIt0(Nx-1+i) = lam2*h(0,i,0,nK2,dx,p,q,a1,a2)*w_pre2;
    end
    FIt0 = dx/2*FIt0;
    % 形成Fct0 (A(j,1))
    Fct0 = zeros(2*Nx-2,1);
    for j =1:Nx-1
        Fct0(j) = lam1*A_K(K,L,q,a2,nK1,nK2,mK1,mK2,dx,dt,1,j,1);
        Fct0(Nx-1+j) = lam2*A_K(K,L,q,a2,nK1,nK2,mK1,mK2,dx,dt,2,j,1);
    end
            
    % 组装F0
    F0t0 = FDt0-dt*FIt0-dt*Fct0;


    %%%%%%%%%%%%%%形成Phi_star_t0

    x = linspace(-L,L,Nx+1);
    w1_star = exp(-nK1*x-mK1*1*dt).*max(K-exp(x),0);
    Phi1_star = w1_star(2:1:Nx).';
    w2_star = exp(-nK2*x-mK2*1*dt).*max(K-exp(x),0);
    Phi2_star = w2_star(2:1:Nx).';
    Phi_star_t0 = [Phi1_star;Phi2_star];


end
function [B,F,Phi_star] = AssembleMatrixMerton(n,dt,dx,Nx,sigma1,sigma2,nM1,nM2,mM1,mM2,A,lam1,lam2,L,K,...
    mu_star1,mu_star2,mu,delt,C1,C2,Phi_pre,Merton1,Merton2)
    % 形成矩阵B,向量F(时间层n层)
   
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
    
    
    
    
    %%%%%%%%%%%%%%形成向量F, F=F0+(A-R)*Phi(n)+B*Phi_star
    %%%%% 形成矩阵A
    
    % 形成矩阵D1,D2
    a1 = ones(1,Nx-1);
    a2 = ones(1,Nx-1);
    for i = 1:Nx-1
        a1(i) = exp(i*dx*(nM2-nM1));
    end
    for i = 1:Nx-1
        a2(i) = exp(i*dx*(nM1-nM2));
    end
    D1 = -A(1,2)*dt*exp(-(nM2-nM1)*L+(mM2-mM1)*(n-1)*dt)*diag(a1);
    D2 = -A(2,1)*dt*exp(-(nM1-nM2)*L+(mM1-mM2)*(n-1)*dt)*diag(a2);
    
    % 组装矩阵A
    A_1 = [-eye(size(B1)) D1; D2 -eye(size(B1))];
    
    % 形成矩阵R
    R = dt*[lam1*Merton1 zeros(size(B1)); zeros(size(B1)) lam2*Merton2];
    
    % 形成向量F0, F0=FD-dt*lam*FI 
    % 形成FD
    w_cur1 = (K-exp(-L))*exp(-mM1*n*dt+nM1*L);
    w_cur2 = (K-exp(-L))*exp(-mM2*n*dt+nM2*L);
    FD = zeros(2*Nx-2,1);
    FD(1) = -beta*sigma1^2/2*w_cur1;     
    FD(Nx) = -beta*sigma2^2/2*w_cur2;
    % 形成FI
    w_pre1 = (K-exp(-L))*exp(-mM1*(n-1)*dt+nM1*L);
    w_pre2 = (K-exp(-L))*exp(-mM2*(n-1)*dt+nM2*L);
    FI = zeros(2*Nx-2,1);
    for i = 1:Nx-1
        FI(i) = f(delt,mu_star1,dx,i,0)*w_pre1;
        FI(Nx-1+i) = f(delt,mu_star2,dx,i,0)*w_pre2;
    end
    FI = dx/2*FI;
    % 形成Fc (A(j,n))
    Fc = zeros(2*Nx-2,1);
    for j =1:Nx-1
        Fc(j) = lam1*A_M(K,L,mu,delt,nM1,nM2,mM1,mM2,dx,dt,1,j,n);
        Fc(Nx-1+j) = lam2*A_M(K,L,mu,delt,nM1,nM2,mM1,mM2,dx,dt,2,j,n);
    end

    % 组装F0
    F0 = FD-dt*FI-dt*Fc;
    % 形成矩阵Phi_star
    x = linspace(-L,L,Nx+1);
    w1_star = exp(-nM1*x-mM1*n*dt).*max(K-exp(x),0);
    Phi1_star = w1_star(2:1:Nx).';
    w2_star = exp(-nM2*x-mM2*n*dt).*max(K-exp(x),0);
    Phi2_star = w2_star(2:1:Nx).';
    Phi_star = [Phi1_star;Phi2_star];
    % 组装F
    %tic
    F = A_1*Phi_pre-R*Phi_pre+F0+B*Phi_star;
    %toc
    
end
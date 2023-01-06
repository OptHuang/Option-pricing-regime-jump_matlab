function [B,F,Phi_star] = AssembleMatrixKou(n,dt,dx,Nx,sigma1,sigma2,nK1,nK2,mK1,mK2,A,lam1,lam2,L,K,...
    p,q,a1,a2,Phi_pre,Kou1,Kou2)
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
    z1 = ones(1,Nx-1);
    z2 = ones(1,Nx-1);
    for i = 1:Nx-1
        z1(i) = exp(i*dx*(nK2-nK1));
        z2(i) = exp(i*dx*(nK1-nK2));
    end
    D1 = -A(1,2)*dt*exp(-(nK2-nK1)*L+(mK2-mK1)*(n-1)*dt)*diag(z1);
    D2 = -A(2,1)*dt*exp(-(nK1-nK2)*L+(mK1-mK2)*(n-1)*dt)*diag(z2);
    
    % 组装矩阵A
    A = [-eye(size(B1)) D1; D2 -eye(size(B1))];
    
    % 形成矩阵R
    R = dt*[lam1*Kou1 zeros(size(B1)); zeros(size(B1)) lam2*Kou2];
    
    % 形成向量F0, F0=FD-dt*lam*FI 
    % 形成FD
    w_cur1 = (K-exp(-L))*exp(-mK1*n*dt+nK1*L);
    w_cur2 = (K-exp(-L))*exp(-mK2*n*dt+nK2*L);
    FD = zeros(2*Nx-2,1);
    FD(1) = -beta*sigma1^2/2*w_cur1;     
    FD(Nx) = -beta*sigma2^2/2*w_cur2;
    % 形成FI
    w_pre1 = (K-exp(-L))*exp(-mK1*(n-1)*dt+nK1*L);
    w_pre2 = (K-exp(-L))*exp(-mK2*(n-1)*dt+nK2*L);
    FI = zeros(2*Nx-2,1);
    for i = 1:Nx-1
        FI(i) = lam1*h(0,i,0,nK1,dx,p,q,a1,a2)*w_pre1;
        FI(Nx-1+i) = lam2*h(0,i,0,nK2,dx,p,q,a1,a2)*w_pre2;
    end
    FI = dx/2*FI;
    % 形成Fc (A(j,n))
    Fc = zeros(2*Nx-2,1);
    for j =1:Nx-1
        Fc(j) = lam1*A_K(K,L,q,a2,nK1,nK2,mK1,mK2,dx,dt,1,j,n);
        Fc(Nx-1+j) = lam2*A_K(K,L,q,a2,nK1,nK2,mK1,mK2,dx,dt,2,j,n);
    end
            
    % 组装F0
    F0 = FD-dt*FI-dt*Fc;
    % 形成矩阵Phi_star
    x = linspace(-L,L,Nx+1);
    w1_star = exp(-nK1*x-mK1*n*dt).*max(K-exp(x),0);
    Phi1_star = w1_star(2:1:Nx).';
    w2_star = exp(-nK2*x-mK2*n*dt).*max(K-exp(x),0);
    Phi2_star = w2_star(2:1:Nx).';
    Phi_star = [Phi1_star;Phi2_star];
    % 组装F
    %tic
    F = A*Phi_pre-R*Phi_pre+F0+B*Phi_star;
    %toc
    
end
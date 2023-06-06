function [B, d1t0, d2t0, F0t0, Phi_star_t0] = AssembleBasicMatMerton(dt,dx,Nx,sigma1,sigma2,nM1,nM2,mM1,mM2,A,lam1,lam2,L,K,...
    mu_star1,mu_star2,mu,delt)

    %beta is the ratio of time-space discretization
    beta = dt/(dx^2);

    %%%%%%%%%%%%%%%%Forming matrix B

    % Forming matrices B1, B2
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

    % Assembling matrix B
    B = [B1 zeros(size(B1)); zeros(size(B1)) B2];

    %%%%%%%%%%%%%%%%Forming diagonal matrix Di elements in vector dit0

    a1 = ones(Nx-1,1);
    a2 = ones(Nx-1,1);
    for i = 1:Nx-1
        a1(i) = exp(i*dx*(nM2-nM1));
        a2(i) = 1/a1(i);
    end
    d1t0 = -A(1,2)*dt*exp(-(nM2-nM1)*L)*a1;
    d2t0 = -A(2,1)*dt*exp(-(nM1-nM2)*L)*a2;

    %%%%%%%%%%%%%%%%Forming F0t0

    % Forming vector F0t0, F0t0=FDt0 -dt*lam*FIt0 -dt*Fct0
    % Forming FDt0
    w_cur1 = (K-exp(-L))*exp(-mM1*1*dt+nM1*L);
    w_cur2 = (K-exp(-L))*exp(-mM2*1*dt+nM2*L);
    FDt0 = zeros(2*Nx-2,1);
    FDt0(1) = -beta*sigma1^2/2*w_cur1;     
    FDt0(Nx) = -beta*sigma2^2/2*w_cur2;
    % Forming FIt0
    w_pre1 = (K-exp(-L))*exp(nM1*L);
    w_pre2 = (K-exp(-L))*exp(nM2*L);
    FIt0 = zeros(2*Nx-2,1);
    for i = 1:Nx-1
        FIt0(i) = lam1*f(delt,mu_star1,dx,i,0)*w_pre1;
        FIt0(Nx-1+i) = lam2*f(delt,mu_star2,dx,i,0)*w_pre2;
    end
    FIt0 = dx/2*FIt0;
    % Forming Fct0 (A(j,1))
    Fct0 = zeros(2*Nx-2,1);
    for j =1:Nx-1
        Fct0(j) = lam1*A_M(K,L,mu,delt,nM1,nM2,mM1,mM2,dx,dt,1,j,1);
        Fct0(Nx-1+j) = lam2*A_M(K,L,mu,delt,nM1,nM2,mM1,mM2,dx,dt,2,j,1);
    end

    % Assembling F0
    F0t0 = FDt0-dt*FIt0-dt*Fct0;

    %%%%%%%%%%%%%%Forming Phi_star_t0

    x = linspace(-L,L,Nx+1);
    w1_star = exp(-nM1*x-mM1*1*dt).*max(K-exp(x),0);
    Phi1_star = w1_star(2:1:Nx).';
    w2_star = exp(-nM2*x-mM2*1*dt).*max(K-exp(x),0);
    Phi2_star = w2_star(2:1:Nx).';
    Phi_star_t0 = [Phi1_star;Phi2_star];

end
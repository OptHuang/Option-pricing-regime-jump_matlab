function ComputeVtrue

    format long
    
    %% Merton
    [T,K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,~,~,~,~,...
    mu,delt,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2,~,~,~,~,~,~,~,~,~] = ParaImput();
    Nt0 = 2000;
    Nx0 = 2560;
    dt0 = T/Nt0;
    dx0 = 2*L/Nx0;
    [Vmertrue1,Vmertrue2] = DataMatrixMerton(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx0,Nt0,dx0,dt0,...
    mu,delt,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2);
    
    save('Vmertrue20002560.mat','Vmertrue1','Vmertrue2')
    
    %% Kou
    [T,K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,~,~,~,~,...
    ~,~,~,~,~,~,~,~,~,~,p,q,a1,a2,~,nK1,nK2,mK1,mK2] = ParaImput();
    Nt0 = 2000;
    Nx0 = 2560;
    dt0 = T/Nt0;
    dx0 = 2*L/Nx0;
    [Vkoutrue1,Vkoutrue2] = DataMatrixKou(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx0,Nt0,dx0,dt0,...
    p,q,a1,a2,nK1,nK2,mK1,mK2);
    
    save('Vkoutrue20002560.mat','Vkoutrue1','Vkoutrue2')
    
end
function Err
    
    %% Merton PCM
    load('../data_truevalue/class=0_d=0.025_lam1=0.25_lam2=0.20_Nx=512_Nt=104858.mat')
    Nxtrue = size(VecMertrue1,2)-1;

    [T,K,sigma1,sigma2,r1,r2,d1,d2,lam1,lam2,gamma,mun,rho0,eps,epsilon,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2,p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = ParaImput(0);
    
    Nx0 = 128;
    Nt0 = round(Nx0^2/2.5);
    
    mult = Nxtrue/Nx0;
    
    U1_star0 = VecMertrue1;
    U2_star0 = VecMertrue2;
    U1_star = zeros(1,Nx0+1);
    U2_star = zeros(1,Nx0+1);
    
    
    for i = 0:Nx0
        U1_star(i+1) = U1_star0(1+i*mult);
        U2_star(i+1) = U2_star0(1+i*mult);
    end

    tic
    
    [Upcm1,Urbf2] = DataMatrixMerton(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx0,Nt0,2*L/Nx0,T/Nt0,...
    mu,delt,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2);

    toc

    Ur1 = Upcm1(end,:);
    Ur2 = Urbf2(end,:);
    
    intF1 = U1_star(1:end)-Ur1(1:end); 
    intF2 = U2_star(1:end)-Ur2(1:end); 
    sqrNorm1 = intF1*intF1'*2*L/(Nx0);
    sqrNorm2 = intF2*intF2'*2*L/(Nx0);
    Norm = sqrt(sqrNorm1+sqrNorm2);
    res = sprintf("%e", Norm)
    
    %% Merton RBF
    load('../data_truevalue/class=0_d=0.025_lam1=0.25_lam2=0.20_Nx=512_Nt=104858.mat')
    Nxtrue = size(VecMertrue1,2)-1;
    
    [T,K,sigma1,sigma2,r1,r2,d1,d2,lam1,lam2,gamma,mun,rho0,eps,epsilon,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2,p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = ParaImput(0);
    
    Nx0 = 128;
    Nt0 = round(Nx0^2/2.5);
    
    mult = Nxtrue/Nx0;
    
    U1_star0 = VecMertrue1;
    U2_star0 = VecMertrue2;
    U1_star = zeros(1,Nx0+1);
    U2_star = zeros(1,Nx0+1);
    
    
    for i = 0:Nx0
        U1_star(i+1) = U1_star0(1+i*mult);
        U2_star(i+1) = U2_star0(1+i*mult);
    end

    tic
    
    [Coefmatrix,A1] = RBFassemblemat(1,0,1/Nt0,Nx0,A,sigma1,sigma2,r1,r2,d1,d2,lam1,lam2,...
    kappaM,kappaK,mu,delt,p,q,a1,a2,L,epsilon);
    [Lmat,Umat] = lu(Coefmatrix);

    [Urbf1,Urbf2] = SolveRBF(Lmat,Umat,A1,L,Nx0,Nt0,K);
    
    toc
    
    Ur1 = Urbf1(end,:);
    Ur2 = Urbf2(end,:);

    intF1 = U1_star(1:end)-Ur1(1:end); 
    intF2 = U2_star(1:end)-Ur2(1:end); 
    sqrNorm1 = intF1*intF1'*2*L/(Nx0);
    sqrNorm2 = intF2*intF2'*2*L/(Nx0);
    Norm = sqrt(sqrNorm1+sqrNorm2);
    res = sprintf("%e", Norm)
    
    %% Kou PCM
    load('../data_truevalue/class=1_d=0.025_lam1=0.25_lam2=0.20_Nx=1024_Nt=419430.mat')
    Nxtrue = size(VecKoutrue1,2)-1;

    [T,K,sigma1,sigma2,r1,r2,d1,d2,lam1,lam2,gamma,mun,rho0,eps,epsilon,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2,p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = ParaImput(1);
    
    Nx0 = 128;
    Nt0 = round(Nx0^2/2.5);
    
    mult = Nxtrue/Nx0;
    
    U1_star0 = VecKoutrue1;
    U2_star0 = VecKoutrue2;
    U1_star = zeros(1,Nx0+1);
    U2_star = zeros(1,Nx0+1);
    
    
    for i = 0:Nx0
        U1_star(i+1) = U1_star0(1+i*mult);
        U2_star(i+1) = U2_star0(1+i*mult);
    end

    tic
    
    [Upcm1,Urbf2] = DataMatrixKou(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx0,Nt0,2*L/Nx0,T/Nt0,...
    p,q,a1,a2,nK1,nK2,mK1,mK2);

    toc
    
    Ur1 = Upcm1(end,:);
    Ur2 = Urbf2(end,:);
                
    intF1 = U1_star(1:end)-Ur1(1:end); 
    intF2 = U2_star(1:end)-Ur2(1:end); 
    sqrNorm1 = intF1*intF1'*2*L/(Nx0);
    sqrNorm2 = intF2*intF2'*2*L/(Nx0);
    Norm = sqrt(sqrNorm1+sqrNorm2);
    res = sprintf("%e", Norm)
    
    %% Kou RBF
    load('../data_truevalue/class=1_d=0.025_lam1=0.25_lam2=0.20_Nx=1024_Nt=419430.mat')
    Nxtrue = size(VecKoutrue1,2)-1;

    [T,K,sigma1,sigma2,r1,r2,d1,d2,lam1,lam2,gamma,mun,rho0,eps,epsilon,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2,p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = ParaImput(1);
    
    Nx0 = 256;
    Nt0 = round(Nx0^2/2.5);
    
    mult = Nxtrue/Nx0;
    
    U1_star0 = VecKoutrue1;
    U2_star0 = VecKoutrue2;
    U1_star = zeros(1,Nx0+1);
    U2_star = zeros(1,Nx0+1);
    
    
    for i = 0:Nx0
        U1_star(i+1) = U1_star0(1+i*mult);
        U2_star(i+1) = U2_star0(1+i*mult);
    end

    tic
    
    [Coefmatrix,A1] = RBFassemblemat(1,1,1/Nt0,Nx0,A,sigma1,sigma2,r1,r2,d1,d2,lam1,lam2,...
    kappaM,kappaK,mu,delt,p,q,a1,a2,L,epsilon);
    [Lmat,Umat] = lu(Coefmatrix);
    
    [Urbf1,Urbf2] = SolveRBF(Lmat,Umat,A1,L,Nx0,Nt0,K);
    
    toc
    
    Ur1 = Urbf1(end,:);
    Ur2 = Urbf2(end,:);
                
    intF1 = U1_star(1:end)-Ur1(1:end); 
    intF2 = U2_star(1:end)-Ur2(1:end); 
    sqrNorm1 = intF1*intF1'*2*L/(Nx0);
    sqrNorm2 = intF2*intF2'*2*L/(Nx0);
    Norm = sqrt(sqrNorm1+sqrNorm2);
    res = sprintf("%e", Norm)
    
end
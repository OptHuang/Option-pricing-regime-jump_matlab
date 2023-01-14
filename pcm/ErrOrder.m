function ErrOrder
    % 求解误差和收敛阶
    
    
    %% Merton型
    clear
    clc
    
    [T,K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,~,~,~,~,...
    mu,delt,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2,~,~,~,~,~,~,~,~,~] = ParaImput();
    Norm1 = zeros(1,5);
    Norm2 = zeros(1,5);
    
    Nt0 = 100000;
    Nx0 = 4;
    % 导入真解
    load('../data_truevalue/Vmertrue100000512_02.mat')
    U1_star0 = Vmertrue1(end,:);
    U2_star0 = Vmertrue2(end,:);

    len = size(U1_star0,2) - 1;


    for i = 1:5
        mult = len/Nx0/2^(i-1);

        [Upcm1,Upcm2] = DataMatrixMerton(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx0*2^(i-1),Nt0,2*L/Nx0/2^(i-1),T/Nt0,...
        mu,delt,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2);
        Ut1 = Upcm1(end,:);
        Ut2 = Upcm2(end,:);
        for j = 1:Nx0*2^(i-1)
            Ut1(j) = U1_star0(1+(j-1)*mult);
            Ut2(j) = U2_star0(1+(j-1)*mult);
        end

        intF1 = Ut1(1:end)-Upcm1(end,1:end); 
        intF2 = Ut2(1:end)-Upcm2(end,1:end); 

        Norm1(i) = sqrt(intF1*intF1'*2*L/(Nx0*2^(i-1)));
        Norm2(i) = sqrt(intF2*intF2'*2*L/(Nx0*2^(i-1)));

    end

    
    N0 = [Nx0 Nx0*2 Nx0*4 Nx0*8 Nx0*16];
    k = 2;
    Standard = -k*log(N0)+1;
    log(Norm1)
    
    plot(log(N0),Standard,'-b',log(N0),log(Norm1),'-ro',log(N0),log(Norm2)-1,'-g*')
    xlabel('$\log(N_x)$','Interpreter','latex','FontSize',16)
    ylabel('log(error)','Interpreter','latex','FontSize',16)
    axis([1.3,4.2,-9.2,-1.7 ])
    legend('order 2','log(error)--regime 1','log(error)--regime 2','Interpreter','latex','FontSize',12)
    title('Convergence order in Example 1','Interpreter','latex','FontSize',22)
    %saveas(gcf,'myfig.jpg')

    
    
    %% Kou型
    clear
    clc
    
    [T,K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,~,~,~,~,...
    ~,~,~,~,~,~,~,~,~,~,p,q,a1,a2,~,nK1,nK2,mK1,mK2] = ParaImput();
    Norm1 = zeros(1,5);
    Norm2 = zeros(1,5);
    
    Nt0 = 100000;
    Nx0 = 8;


    % 导入真解
    load('../data_truevalue/Vkoutrue100000512_02.mat')
    U1_star0 = Vkoutrue1(end,:);
    U2_star0 = Vkoutrue2(end,:);

    len = size(U1_star0,2) - 1;

    for i = 1:5
        mult = len/Nx0/2^(i-1);

        [Upcm1,Upcm2] = DataMatrixKou(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx0*2^(i-1),Nt0,2*L/Nx0/2^(i-1),T/Nt0,...
        p,q,a1,a2,nK1,nK2,mK1,mK2);
        Ut1 = Upcm1(end,:);
        Ut2 = Upcm2(end,:);
        for j = 1:Nx0*2^(i-1)
            Ut1(j) = U1_star0(1+(j-1)*mult);
            Ut2(j) = U2_star0(1+(j-1)*mult);
        end

        intF1 = Ut1(1:end)-Upcm1(end,1:end); 
        intF2 = Ut2(1:end)-Upcm2(end,1:end); 

        Norm1(i) = sqrt(intF1*intF1'*2*L/(Nx0*2^(i-1)));
        Norm2(i) = sqrt(intF2*intF2'*2*L/(Nx0*2^(i-1)));
        
    end
    
    N0 = [Nx0 Nx0*2 Nx0*4 Nx0*8 Nx0*16];
    k = 2;
    Standard = -k*log(N0)+2;
    log(Norm1)
    
    plot(log(N0),Standard,'-b',log(N0),log(Norm1),'-ro',log(N0),log(Norm2)-1,'-g*')
    xlabel('$\log(N_x)$','Interpreter','latex','FontSize',16)
    ylabel('log(error)','Interpreter','latex','FontSize',16)
    axis([2,5,-10,-2 ])
    legend('order 2','log(error)--regime 1','log(error)--regime 2','Interpreter','latex','FontSize',12)
    title('Convergence order in Example 2','Interpreter','latex','FontSize',22)
    %saveas(gcf,'myfig.jpg')

end
function CompareFreeBoundary
    % 比较不同几种期权的自由边界
    
    %% Merton
    [T,K,sigma1,sigma2,r1,r2,lam1,lam2,gamma,mun,rho0,eps,epsilon,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2,p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = ParaImput();
    
    t = linspace(0,1,Nt+1);
    x = linspace(-L,L,Nx+1);
    s = exp(x);
    
    % 计算两种波动率下的永久美式看跌跳扩散自由边界 B_bar1
    max_iter = 100;
    eps = 1e-9;
    beta1 = -Newton_perMerton(max_iter,eps,sigma1,r1,0,lam1,kappaM,mu,delt);
    B_bar1 = K*beta1/(1+beta1);

     % 计算两种波动率下的美式看跌跳扩散的自由边界曲线 B_star1,B_star2
    A_star = 0*A;
    [mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2] = Mertonpara(sigma1,sigma2,r1,r2,0,0,lam1,lam2,A_star);
    [Vstar1,Vstar2] = DataMatrixMerton(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A_star,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2);
    [~,~,fb1,fb2] = FreeBoundry(Vstar1,Vstar2,K,L,Nx,Nt);
    B_star1 = K-flipud(fb1);
    B_star2 = K-flipud(fb2);
    
    % 计算两种波动率下的美式看跌制度转换跳扩散 B_1,B_2
    [mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2] = Mertonpara(sigma1,sigma2,r1,r2,0,0,lam1,lam2,A);
    [Vtrue1,Vtrue2] = DataMatrixMerton(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2);
    [~,~,fb_1,fb_2] = FreeBoundry(Vtrue1,Vtrue2,K,L,Nx,Nt);
    B_1 = K-flipud(fb_1);
    B_2 = K-flipud(fb_2);
    
    % 绘制图像
    plot(t,B_2,'b-')
    hold on
    plot(t,B_1,'g-')
    hold on
    plot(t,B_star1,'r-')
    hold on
    plot(t,B_bar1*ones(size(t)),'Color','#0072BD','linewidth',3)
    hold on

    
    ylabel('$S$','interpreter','latex','FontSize',16)
    xlabel('$t$','interpreter','latex','FontSize',16)
    ylim([0 1])
    xlim([0 1])
    set(legend('$B_{2}(t)$','$B_{1}(t)$','$\tilde{B}_{1}(t)$','$\bar{B}_{1}$',Location='northwest'),...
        'interpreter','latex','FontSize',12)
    title('Optimal exercise boundaries in Example 1','Interpreter','latex','FontSize',22)

    
    %% Kou
    [T,K,sigma1,sigma2,r1,r2,lam1,lam2,gamma,mun,rho0,eps,epsilon,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2,p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = ParaImput();

    t = linspace(0,1,Nt+1);
    x = linspace(-L,L,Nx+1);
    s = exp(x);
    
    % 计算两种波动率下的永久美式看跌跳扩散自由边界 B_bar1
    omg1 = r1-0-lam1*kappaK-0.5*sigma1^2;
    para1 = [-0.5*sigma1^2,(0.5*sigma1^2*(a1-a2)-omg1),(0.5*sigma1*a1*a2+omg1*(a1-a2)+(r1+lam1)),...
        (omg1*a1*a2-(r1+lam1)*(a1-a2)+lam1*(p*a1-q*a2)),-r1*a1*a2];
    x1 = roots(para1);  % 求四次多项式的根
    x1 = sort(x1);
    beta41 = -x1(1);
    beta31 = -x1(2);
    B_bar1 = K*(a2+1)/a2*beta31/(1+beta31)*beta41/(1+beta41);

     % 计算两种波动率下的美式看跌跳扩散的自由边界曲线 B_star1,B_star2
    A_star = 0*A;
    [p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = Koupara(sigma1,sigma2,r1,r2,0,0,lam1,lam2,A_star);
    [Vstar1,Vstar2] = DataMatrixKou(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A_star,x0,t0,Nx,Nt,dx,dt,...
    p,q,a1,a2,nK1,nK2,mK1,mK2);
    [~,~,fb1,fb2] = FreeBoundry(Vstar1,Vstar2,K,L,Nx,Nt);
    B_star1 = K-flipud(fb1);
    B_star2 = K-flipud(fb2);
    
    % 计算两种波动率下的美式看跌制度转换跳扩散 B_1,B_2
    [p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = Koupara(sigma1,sigma2,r1,r2,0,0,lam1,lam2,A);
    [Vtrue1,Vtrue2] = DataMatrixKou(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx,Nt,dx,dt,...
    p,q,a1,a2,nK1,nK2,mK1,mK2);
    [~,~,fb_1,fb_2] = FreeBoundry(Vtrue1,Vtrue2,K,L,Nx,Nt);
    B_1 = K-flipud(fb_1);
    B_2 = K-flipud(fb_2);
    
    % 绘制图像g
    plot(t,B_2,'b-')
    hold on
    plot(t,B_1,'g-')
    hold on
    plot(t,B_star1,'r-')
    hold on
    plot(t,B_bar1*ones(size(t)),'Color','#0072BD','linewidth',3)
    hold on

 
    ylabel('$S$','interpreter','latex','FontSize',16)
    xlabel('$t$','interpreter','latex','FontSize',16)
    ylim([0 1])
    xlim([0 1])
    set(legend('$B_{2}(t)$','$B_{1}(t)$','$\tilde{B}_{1}(t)$','$\bar{B}_{1}$',Location='northwest'),...
        'interpreter','latex','FontSize',12)
    title('Optimal exercise boundaries in Example 2','Interpreter','latex','FontSize',22)

end
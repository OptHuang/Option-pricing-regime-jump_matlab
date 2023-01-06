function Pic3DPCM
    % 绘制三维曲面(PCM)
    
    
    %% Merton
    format long
    tic
    
    [T,K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2,~,~,~,~,~,~,~,~,~] = ParaImput();
    
    [Vtrue1,Vtrue2] = DataMatrixMerton(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2);
    [~,~,fb1,fb2] = FreeBoundry(Vtrue1,Vtrue2,K,L,Nx,Nt);
    [S,M,t,~] = MeshGeneration(T,L,Nx,Nt);
    
%   % 现实中矩阵元素对应坐标需要将矩阵以空间方向为轴翻转一下
    Vtrue1 = flipud(Vtrue1);
    Vtrue2 = flipud(Vtrue2);
    fb1 = flipud(fb1);
    fb2 = flipud(fb2);
    
    subplot(1,2,1)
    mesh(S,M,Vtrue1)
    hold on
    plot3(K-fb1,t,fb1,'k-','LineWidth',3)
    hold on
    plot3(K-fb1,t,zeros(size(t)),'k-','LineWidth',3)
    axis([0 10 0 1 0 1.2])
    xlabel('S')
    ylabel('T')
    zlabel('P')
    title('Example1: option value under regime 1')
    
    subplot(1,2,2)
    mesh(S,M,Vtrue2)
    hold on
    plot3(K-fb2,t,fb2,'k-','LineWidth',3)
    hold on
    plot3(K-fb1,t,zeros(size(t)),'k-','LineWidth',3)
    axis([0 10 0 1 0 1.2])
    xlabel('S')
    ylabel('T')
    zlabel('P')
    title('Example1: option value under regime 2')
    
%     subplot(1,3,3)
%     x = linspace(-L,L,Nx+1);
%     s = exp(x);
%     plot(s,Utrue1(1,:))
%     axis([0 25 0 9])
    
    toc
    
    %% Kou
    format long
    tic
    
    [T,K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx,Nt,dx,dt,...
    ~,~,~,~,~,~,~,~,~,~,p,q,a1,a2,~,nK1,nK2,mK1,mK2] = ParaImput();
    
    [Vtrue1,Vtrue2] = DataMatrixKou(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx,Nt,dx,dt,...
    p,q,a1,a2,nK1,nK2,mK1,mK2);
    [~,~,fb1,fb2] = FreeBoundry(Vtrue1,Vtrue2,K,L,Nx,Nt);
    [S,M,t,~] = MeshGeneration(T,L,Nx,Nt);
    
%   % 现实中矩阵元素对应坐标需要将矩阵以空间方向为轴翻转一下
    Vtrue1 = flipud(Vtrue1);
    Vtrue2 = flipud(Vtrue2);
    fb1 = flipud(fb1);
    fb2 = flipud(fb2);
    
    subplot(1,2,1)
    mesh(S,M,Vtrue1)
    hold on
    plot3(K-fb1,t,fb1,'k-','LineWidth',3)
    hold on
    plot3(K-fb1,t,zeros(size(t)),'k-','LineWidth',3)
    axis([0 10 0 1 0 1.2])
    xlabel('S')
    ylabel('T')
    zlabel('P')
    title('Example1: option value under regime 1')
    
    subplot(1,2,2)
    mesh(S,M,Vtrue2)
    hold on
    plot3(K-fb2,t,fb2,'k-','LineWidth',3)
    hold on
    plot3(K-fb1,t,zeros(size(t)),'k-','LineWidth',3)
    axis([0 10 0 1 0 1.2])
    xlabel('S')
    ylabel('T')
    zlabel('P')
    title('Example1: option value under regime 2')
    
%     subplot(1,3,3)
%     x = linspace(-L,L,Nx+1);
%     s = exp(x);
%     plot(s,Utrue1(1,:))
%     axis([0 25 0 9])
    
    toc
    
end
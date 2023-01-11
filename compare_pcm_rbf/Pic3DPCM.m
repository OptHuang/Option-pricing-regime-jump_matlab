function Pic3DPCM
    % 绘制三维曲面(PCM)
    
    
    %% Merton
    format long
    tic
    
    [T,K,sigma1,sigma2,r1,r2,lam1,lam2,gamma,mun,rho0,eps,epsilon,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2,p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = ParaImput();
    
    [Vtrue1,Vtrue2] = DataMatrixMerton(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2);
    [~,~,fb1,fb2] = FreeBoundry(Vtrue1,Vtrue2,K,L,Nx,Nt);
    [S,M,t,~] = MeshGeneration(T,L,Nx,Nt);
    
    % 现实中矩阵元素对应坐标需要将矩阵以空间方向为轴翻转一下
    Vtrue1 = flipud(Vtrue1);
    Vtrue2 = flipud(Vtrue2);
    fb1 = flipud(fb1);
    fb2 = flipud(fb2);
    
    % 补充S=0
    V1true = zeros(Nt+1,Nx+2);
    V2true = zeros(Nt+1,Nx+2);
    V1true(:,1) = K;
    V2true(:,1) = K;
    V1true(:,2:end) = Vtrue1;
    V2true(:,2:end) = Vtrue2;

%     subplot(1,1,1)
%     mesh(S,M,V1true)
%     hold on
%     plot3(K-fb1,t,fb1,'k-','LineWidth',3)
%     %hold on
%     %plot3(K-fb1,t,zeros(size(t)),'k-','LineWidth',3)
%     axis([0 7.5 0 1 0 1])
%     xlabel('$S$','Interpreter','latex','FontSize',16)
%     ylabel('$t$','Interpreter','latex','FontSize',16)
%     zlabel('$P$','Interpreter','latex','FontSize',16)
%     title('Option price under regime 1 in Example 1','Interpreter','latex','FontSize',22)
%     set(gcf,'position',[10,5,560,420])
    
    
    subplot(1,1,1)
    mesh(S,M,V2true)
    hold on
    plot3(K-fb2,t,fb2,'k-','LineWidth',3)
    %hold on
    %plot3(K-fb1,t,zeros(size(t)),'k-','LineWidth',3)
    axis([0 7.5 0 1 0 1])
    xlabel('$S$','Interpreter','latex','FontSize',16)
    ylabel('$t$','Interpreter','latex','FontSize',16)
    zlabel('$P$','Interpreter','latex','FontSize',16)
    title('Option price under regime 2 in Example 1','Interpreter','latex','FontSize',22)
    set(gcf,'position',[10,5,560,420])
  
    
    toc
    
    %% Kou
    format long
    tic
    
    [T,K,sigma1,sigma2,r1,r2,lam1,lam2,gamma,mun,rho0,eps,epsilon,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2,p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = ParaImput();
    
    [Vtrue1,Vtrue2] = DataMatrixKou(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx,Nt,dx,dt,...
    p,q,a1,a2,nK1,nK2,mK1,mK2);
    [~,~,fb1,fb2] = FreeBoundry(Vtrue1,Vtrue2,K,L,Nx,Nt);
    [S,M,t,~] = MeshGeneration(T,L,Nx,Nt);
    
    % 现实中矩阵元素对应坐标需要将矩阵以空间方向为轴翻转一下
    Vtrue1 = flipud(Vtrue1);
    Vtrue2 = flipud(Vtrue2);
    fb1 = flipud(fb1);
    fb2 = flipud(fb2);

    % 补充S=0
    V1true = zeros(Nt+1,Nx+2);
    V2true = zeros(Nt+1,Nx+2);
    V1true(:,1) = K;
    V2true(:,1) = K;
    V1true(:,2:end) = Vtrue1;
    V2true(:,2:end) = Vtrue2;
    
    subplot(1,1,1)
    mesh(S,M,V1true)
    hold on
    plot3(K-fb1,t,fb1,'k-','LineWidth',3)
    %hold on
    %plot3(K-fb1,t,zeros(size(t)),'k-','LineWidth',3)
    axis([0 7.5 0 1 0 1])
    xlabel('$S$','Interpreter','latex','FontSize',16)
    ylabel('$t$','Interpreter','latex','FontSize',16)
    zlabel('$P$','Interpreter','latex','FontSize',16)
    title('Option price under regime 1 in Example 2','Interpreter','latex','FontSize',22)
    set(gcf,'position',[10,5,560,420])

    
%     subplot(1,1,1)
%     mesh(S,M,V2true)
%     hold on
%     plot3(K-fb2,t,fb2,'k-','LineWidth',3)
%     %hold on
%     %plot3(K-fb1,t,zeros(size(t)),'k-','LineWidth',3)
%     axis([0 7.5 0 1 0 1])
%     xlabel('$S$','Interpreter','latex','FontSize',16)
%     ylabel('$t$','Interpreter','latex','FontSize',16)
%     zlabel('$P$','Interpreter','latex','FontSize',16)
%     title('Option price under regime 2 in Example 2','Interpreter','latex','FontSize',22)
%     set(gcf,'position',[10,5,560,420])
    
    toc
    
end
function contract_pcmrbf_samefig
    % 在一张图上比较pcm和rbf算出来的解, 看是否吻合
    
    %% Merton
    [T,K,sigma1,sigma2,r1,r2,lam1,lam2,gamma,mun,rho0,eps,epsilon,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2,p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = ParaImput();

    x = linspace(-L,L,Nx+1); 
    s = zeros(1,Nx+2);
    s(1) = 0;
    s(2:end) = exp(x);

    [Vpcmtrue1,Vpcmtrue2] = DataMatrixMerton(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2);
    Vpcm1 = zeros(1,Nx+2);
    Vpcm2 = zeros(1,Nx+2);
    Vpcm1(1) = 1;
    Vpcm2(1) = 1;
    Vpcm1(2:end) = Vpcmtrue1(end,:);
    Vpcm2(2:end) = Vpcmtrue2(end,:);
    
    [Coefmatrix,A1] = RBFassemblemat(1,0,dt,Nx,A,sigma1,sigma2,r1,r2,lam1,lam2,...
    kappaM,kappaK,mu,delt,p,q,a1,a2,L,epsilon);
    [Lmat,Umat] = lu(Coefmatrix);
    [Vrbftrue1,Vrbftrue2] = SolveRBF(Lmat,Umat,A1,L,Nx,Nt,K);
    Vrbf1 = zeros(1,Nx+2);
    Vrbf2 = zeros(1,Nx+2);
    Vrbf1(1) = 1;
    Vrbf2(1) = 1;
    Vrbf1(2:end) = Vrbftrue1(end,:);
    Vrbf2(2:end) = Vrbftrue2(end,:);

    Vbound = max(K-s,0);
    
    % 隔8个点选一个点 (因为圆圈太密集)
    vpcm1 = zeros(1,Nx/8+2);
    vpcm2 = zeros(1,Nx/8+2);
    s1 = zeros(1,Nx/8+2);
    vpcm1(1,1) = Vpcm1(1,1);
    vpcm2(1,1) = Vpcm2(1,1);
    s1(1,1) = s(1,1);
    vpcm1(1,2) = Vpcm1(1,2);
    vpcm2(1,2) = Vpcm2(1,2);
    s1(1,2) = s(1,2);
    for i = 1:Nx/8
        vpcm1(1,i+2) = Vpcm1(1,8*i+2);
        vpcm2(1,i+2) = Vpcm2(1,8*i+2);
        s1(1,i+2) = s(1,8*i+2);
    end

    % 子图采点
    smin = 1.2;
    smax = 2.4;
    x_sub_left = max(find(s<smin));
    x_sub_right = min(find(s>smax));
    x_sub = s(x_sub_left:x_sub_right);
    yrbf1_sub = Vrbf1(x_sub_left:x_sub_right);
    yrbf2_sub = Vrbf2(x_sub_left:x_sub_right);
    x1_sub_left = max(find(s1<smin));
    x1_sub_right = min(find(s1>smax));
    x1_sub = s1(x1_sub_left:x1_sub_right);
    ypcm1_sub = vpcm1(x1_sub_left:x1_sub_right);
    ypcm2_sub = vpcm2(x1_sub_left:x1_sub_right);
    
    % 主图像
    plot(s1,vpcm1,'ro')
    hold on
    plot(s1,vpcm2,'bo')
    hold on
    plot(s,Vrbf1,'m-')
    hold on
    plot(s,Vrbf2,'c-')
    hold on
    plot(s,Vbound,'g--','linewidth',1)
    
    xlim([0,4])
    xlabel('$S$','Interpreter','latex')
    ylabel('$P$','Interpreter','latex')
    title('Option price at $t=0$ in Example 1','Interpreter','latex')
    legend({'FDPCM--$P_{1}(S,0)$','FDPCM--$P_{2}(S,0)$','RBF--$P_{1}(S,0)$','RBF--$P_{2}(S,0)$',...
        '$\max\{K-S,0\}$'},'Interpreter','latex')

    % 红色矩形框
    rectangle('Position',[smin 0.01 smax-smin 0.17], 'EdgeColor', 'r')
    hold on

    % 子图像
    H = axes('Position',[0.35,0.42,0.28,0.25]); % 生成子图
    plot(x_sub,yrbf1_sub,'m-',x_sub,yrbf2_sub,'c-',x1_sub,ypcm1_sub,'ro',x1_sub,ypcm2_sub,'bo');                % 绘制局部曲线图
    xlim([min(x_sub),max(x_sub)]);    % 设置坐标轴范围    
    set(H, 'XTick',[], 'YTick', []);

    
    
    
    %% Kou
    [T,K,sigma1,sigma2,r1,r2,lam1,lam2,gamma,mun,rho0,eps,epsilon,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2,p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = ParaImput();

    x = linspace(-L,L,Nx+1); 
    s = zeros(1,Nx+2);
    s(1) = 0;
    s(2:end) = exp(x);

    [Vpcmtrue1,Vpcmtrue2] = DataMatrixKou(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx,Nt,dx,dt,...
    p,q,a1,a2,nK1,nK2,mK1,mK2);
    Vpcm1 = zeros(1,Nx+2);
    Vpcm2 = zeros(1,Nx+2);
    Vpcm1(1) = 1;
    Vpcm2(1) = 1;
    Vpcm1(2:end) = Vpcmtrue1(end,:);
    Vpcm2(2:end) = Vpcmtrue2(end,:);
    
    [Coefmatrix,A1] = RBFassemblemat(1,1,dt,Nx,A,sigma1,sigma2,r1,r2,lam1,lam2,...
    kappaM,kappaK,mu,delt,p,q,a1,a2,L,epsilon);
    [Lmat,Umat] = lu(Coefmatrix);
    [Vrbftrue1,Vrbftrue2] = SolveRBF(Lmat,Umat,A1,L,Nx,Nt,K);
    Vrbf1 = zeros(1,Nx+2);
    Vrbf2 = zeros(1,Nx+2);
    Vrbf1(1) = 1;
    Vrbf2(1) = 1;
    Vrbf1(2:end) = Vrbftrue1(end,:);
    Vrbf2(2:end) = Vrbftrue2(end,:);
    
    Vbound = max(K-s,0);
    
    % 隔8个点选一个点 (因为圆圈太密集)
    vpcm1 = zeros(1,Nx/8+2);
    vpcm2 = zeros(1,Nx/8+2);
    s1 = zeros(1,Nx/8+2);
    vpcm1(1,1) = Vpcm1(1,1);
    vpcm2(1,1) = Vpcm2(1,1);
    s1(1,1) = s(1,1);
    vpcm1(1,2) = Vpcm1(1,2);
    vpcm2(1,2) = Vpcm2(1,2);
    s1(1,2) = s(1,2);
    for i = 1:Nx/8
        vpcm1(1,i+2) = Vpcm1(1,8*i+2);
        vpcm2(1,i+2) = Vpcm2(1,8*i+2);
        s1(1,i+2) = s(1,8*i+2);
    end

    % 子图采点
    smin = 1.2;
    smax = 2.4;
    x_sub_left = max(find(s<smin));
    x_sub_right = min(find(s>smax));
    x_sub = s(x_sub_left:x_sub_right);
    yrbf1_sub = Vrbf1(x_sub_left:x_sub_right);
    yrbf2_sub = Vrbf2(x_sub_left:x_sub_right);
    x1_sub_left = max(find(s1<smin));
    x1_sub_right = min(find(s1>smax));
    x1_sub = s1(x1_sub_left:x1_sub_right);
    ypcm1_sub = vpcm1(x1_sub_left:x1_sub_right);
    ypcm2_sub = vpcm2(x1_sub_left:x1_sub_right);

    % 主图像
    plot(s1,vpcm1,'ro')
    hold on
    plot(s1,vpcm2,'bo')
    hold on
    plot(s,Vrbf1,'m-')
    hold on
    plot(s,Vrbf2,'c-')
    hold on
    plot(s,Vbound,'g--','linewidth',1)

    xlim([0,4])
    xlabel('$S$','Interpreter','latex')
    ylabel('$P$','Interpreter','latex')
    title('Option price at $t=0$ in Example 2','Interpreter','latex')
    legend({'FDPCM--$P_{1}(S,0)$','FDPCM--$P_{2}(S,0)$','RBF--$P_{1}(S,0)$','RBF--$P_{2}(S,0)$',...
        '$\max\{K-S,0\}$'},'Interpreter','latex')

    % 红色矩形框
    rectangle('Position',[smin 0.01 smax-smin 0.19], 'EdgeColor', 'r')
    hold on

    % 子图像
    H = axes('Position',[0.35,0.42,0.28,0.25]); % 生成子图
    plot(x_sub,yrbf1_sub,'m-',x_sub,yrbf2_sub,'c-',x1_sub,ypcm1_sub,'ro',x1_sub,ypcm2_sub,'bo');                % 绘制局部曲线图
    xlim([min(x_sub),max(x_sub)]);    % 设置坐标轴范围    
    set(H, 'XTick',[], 'YTick', []);



end
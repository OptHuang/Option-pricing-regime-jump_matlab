function Plot2D_compare(Nx,Nt,problem_paras,options)

    currentFolder = fileparts(mfilename('fullpath'));
    testSubFolder = fullfile(currentFolder, '..', 'methods');
    addpath(testSubFolder);

    Vpcm1 = zeros(1,Nx+2);
    Vpcm2 = zeros(1,Nx+2);
    Vrbf1 = zeros(1,Nx+2);
    Vrbf2 = zeros(1,Nx+2);
    Vpcm1(1) = problem_paras.K;
    Vpcm2(1) = problem_paras.K;
    Vrbf1(1) = problem_paras.K;
    Vrbf2(1) = problem_paras.K;

    if nargin < 4
        [last_row_pcm1,last_row_pcm2,~,~,~,~,L] = FDPCM(Nx,Nt,problem_paras);
        [last_row_rbf1,last_row_rbf2,~,~,~,~,~] = RBCM(Nx,Nt,problem_paras);
    else
        [last_row_pcm1,last_row_pcm2,~,~,~,~,L] = FDPCM(Nx,Nt,problem_paras,options);
        [last_row_rbf1,last_row_rbf2,~,~,~,~,~] = RBCM(Nx,Nt,problem_paras,options);
    end

    x = linspace(-L,L,Nx+1); 
    s = zeros(1,Nx+2);
    s(1) = 0;
    s(2:end) = exp(x);

    Vpcm1(2:end) = last_row_pcm1;
    Vpcm2(2:end) = last_row_pcm2;
    Vrbf1(2:end) = last_row_rbf1;
    Vrbf2(2:end) = last_row_rbf2;

    Vbound = max(problem_paras.K-s,0);

    % 隔8个点选一个点 (因为圆圈太密集)
    vpcm1 = zeros(1,floor(Nx/8)+2);
    vpcm2 = zeros(1,floor(Nx/8)+2);
    s1 = zeros(1,floor(Nx/8)+2);
    vpcm1(1,1) = Vpcm1(1,1);
    vpcm2(1,1) = Vpcm2(1,1);
    s1(1,1) = s(1,1);
    vpcm1(1,2) = Vpcm1(1,2);
    vpcm2(1,2) = Vpcm2(1,2);
    s1(1,2) = s(1,2);
    for i = 1:floor(Nx/8)
        vpcm1(1,i+2) = Vpcm1(1,8*i+2);
        vpcm2(1,i+2) = Vpcm2(1,8*i+2);
        s1(1,i+2) = s(1,8*i+2);
    end

    % 子图采点
    smin = 1.2;
    smax = 2.4;
    x_sub_left = find(s<smin, 1, 'last' );
    x_sub_right = find(s>smax, 1 );
    x_sub = s(x_sub_left:x_sub_right);
    yrbf1_sub = Vrbf1(x_sub_left:x_sub_right);
    yrbf2_sub = Vrbf2(x_sub_left:x_sub_right);
    x1_sub_left = find(s1<smin, 1, 'last' );
    x1_sub_right = find(s1>smax, 1 );
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
    xlabel('$S$','Interpreter','latex','FontSize',16)
    ylabel('$P$','Interpreter','latex','FontSize',16)

    switch problem_paras.type
        case {"Merton"}
            title('Option price at $t=0$ in Example 1','Interpreter','latex','FontSize',22)
        case {"Kou"}
            title('Option price at $t=0$ in Example 2','Interpreter','latex','FontSize',22)
    end
    
    legend({'FDPCM--$P_{1}(S,0)$','FDPCM--$P_{2}(S,0)$','RBCM--$P_{1}(S,0)$','RBCM--$P_{2}(S,0)$',...
        '$\max\{K-S,0\}$'},'Interpreter','latex','FontSize',12)

    % 红色矩形框
    rectangle('Position',[smin 0.01 smax-smin-0.2 0.17], 'EdgeColor', 'r')
    hold on

    % 子图像
    H = axes('Position',[0.35,0.42,0.26,0.25]); % 生成子图
    plot(x_sub,yrbf1_sub,'m-',x_sub,yrbf2_sub,'c-',x1_sub,ypcm1_sub,'ro',x1_sub,ypcm2_sub,'bo');             
    
    % 绘制局部曲线图
    xlim([min(x_sub),max(x_sub)-0.2]);    % 设置坐标轴范围    
    set(H, 'XTick',[], 'YTick', []);

end
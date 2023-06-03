function Plot_Freeboundaries(Nx,Nt,problem_paras,options)

    currentFolder = fileparts(mfilename('fullpath'));
    testSubFolder = fullfile(currentFolder, '..', 'methods');
    addpath(testSubFolder);
    FDPCM_subfuncFolder = fullfile(currentFolder, '..', 'methods/FDPCM_subfunc');
    addpath(FDPCM_subfuncFolder)

    t = linspace(0,1,Nt+1);

    % 计算两种波动率下的永久美式看跌跳扩散自由边界 B_bar1
    switch problem_paras.type
        case {"Merton"}
            max_iter = 100;
            eps = 1e-9;
            [mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2] = ...
                Mertonpara(problem_paras.sig1,problem_paras.sig2,problem_paras.r1,...
                problem_paras.r2,problem_paras.d1,problem_paras.d2,problem_paras.lam1,...
                problem_paras.lam2,problem_paras.A);
            beta1 = -Newton_perMerton(max_iter,eps,problem_paras.sig1,...
                problem_paras.r1,problem_paras.d1,problem_paras.lam1,kappaM,mu,delt);
            B_bar1 = problem_paras.K*beta1/(1+beta1);
        case {"Kou"}
            [p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = Koupara(problem_paras.sig1,...
                problem_paras.sig2,problem_paras.r1,problem_paras.r2,...
                problem_paras.d1,problem_paras.d2,problem_paras.lam1,...
                problem_paras.lam2,problem_paras.A);
            omg1 = problem_paras.r1-problem_paras.d1-...
                problem_paras.lam1*kappaK-0.5*problem_paras.sig1^2;
            para1 = [-0.5*problem_paras.sig1^2,...
                (0.5*problem_paras.sig1^2*(a1-a2)-omg1),...
                (0.5*problem_paras.sig1*a1*a2+omg1*(a1-a2)+...
                (problem_paras.r1+problem_paras.lam1)),(omg1*a1*a2-...
                (problem_paras.r1+problem_paras.lam1)*...
                (a1-a2)+problem_paras.lam1*(p*a1-q*a2)),...
                -problem_paras.r1*a1*a2];
            x1 = roots(para1);  % 求四次多项式的根
            x1 = sort(x1);
            beta41 = -x1(1);
            beta31 = -x1(2);
            B_bar1 = problem_paras.K*(a2+1)/a2*beta31/(1+beta31)*beta41/(1+beta41);
    end
    
    if nargin < 4
        % 计算两种波动率下的美式看跌制度转换跳扩散 B_1,B_2
        [~,~,~,~,fb_1,fb_2,~] = FDPCM(Nx,Nt,problem_paras);
        B_1 = problem_paras.K-flipud(fb_1);
        B_2 = problem_paras.K-flipud(fb_2);
    
        % 计算两种波动率下的美式看跌跳扩散的自由边界曲线 B_star1
        problem_paras.A = 0*problem_paras.A;
        [~,~,~,~,fb1,~,~] = FDPCM(Nx,Nt,problem_paras);
        B_star1 = problem_paras.K-flipud(fb1);
    else
        % 计算两种波动率下的美式看跌制度转换跳扩散 B_1,B_2
        [~,~,~,~,fb_1,fb_2,~] = FDPCM(Nx,Nt,problem_paras,options);
        B_1 = problem_paras.K-flipud(fb_1);
        B_2 = problem_paras.K-flipud(fb_2);
    
        % 计算两种波动率下的美式看跌跳扩散的自由边界曲线 B_star1
        problem_paras.A = 0*problem_paras.A;
        [~,~,~,~,fb1,~,~] = FDPCM(Nx,Nt,problem_paras,options);
        B_star1 = problem_paras.K-flipud(fb1);
    end
    
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

    switch problem_paras.type
        case {"Merton"}
            title('Optimal exercise boundaries in Example 1','Interpreter','latex','FontSize',22)
        case {"Kou"}
            title('Optimal exercise boundaries in Example 2','Interpreter','latex','FontSize',22)
    end

end
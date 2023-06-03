function Phi_end = PCMPro(eps,nu,mu,vrho,M,F,Phi_start)
    % pcm算法核心程序

    % Nx是空间方向剖分段数
    % eps是终止迭代误差
    % nu、mu、vrho是给定的(0,1)之间的系数
    % M是给定的(2Nx-2)*(2Nx-2)的系数矩阵
    % F是给定的(2Nx-2)*1的列向量
    % Phi_start是该步迭代的初始值
    
    len = size(F,1);
    Phi = Phi_start;
    
    beta = 1;
    k = 0;
    F_u = M'*Phi+F;
    tol = norm(Phi-max(Phi-F_u,zeros(len,1)),inf);
    while (tol>eps&&k<10000)
        % 投影
        Phi_mid = Phi;
        F_u_mid = F_u;
        Phi = max(Phi_mid-beta.*F_u_mid,0);
        % 可行方向
        F_u = M'*Phi+F;
        du = Phi_mid-Phi;
        dF = beta.*(F_u_mid-F_u);
        rho = norm(dF)/norm(du);
        while (rho>nu)
            beta = 2/3*beta*min(1,1/rho);
            Phi = max(Phi_mid-beta.*F_u_mid,zeros(len,1));
            F_u = M'*Phi+F;
            du = Phi_mid-Phi;
            dF = beta*(F_u_mid-F_u);
            rho = norm(dF)/norm(du);
        end
        duF = du+dF;
        r1 = dot(du,du);
        r2 = dot(duF,duF);
        alpha = r1/r2;
        Phi = Phi_mid-alpha*vrho*duF;
        F_u = M'*Phi+F;
        tol = norm(Phi-max(Phi-F_u,zeros(len,1)),inf);
        if (rho<mu)
           beta = beta*vrho; 
        end
        k = k+1;
    end     
    Phi_end = Phi;
    
end
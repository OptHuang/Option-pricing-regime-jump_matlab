function Phi_end = PCMPro(Nx,eps,gamma,mu,rho0,M,F,Phi_start)
    % pcm算法核心程序

    % Nx是空间方向剖分段数
    % eps是终止迭代误差
    % gamma、mu、rho0是给定的(0,1)之间的系数
    % M是给定的(2Nx-2)*(2Nx-2)的系数矩阵
    % F是给定的(2Nx-2)*1的列向量
    % Phi_start是该步迭代的初始值
    
    Phi = Phi_start;
    
    beta = 1;
    k = 0;
    F_u = M'*Phi+F;
    tol = norm(Phi-max(Phi-F_u,zeros(2*Nx-2,1)),inf);
    while (tol>eps&&k<5000)
        % 投影
        Phi_mid = Phi;
        F_u_mid = F_u;
        Phi = max(Phi_mid-beta.*F_u_mid,0);
        % 可行方向
        F_u = M'*Phi+F;
        du = Phi_mid-Phi;
        dF = beta.*(F_u_mid-F_u);
        rho = norm(dF)/norm(du);
        while (rho>gamma)
            beta = 2/3*beta*min(1,1/rho);
            Phi = max(Phi_mid-beta.*F_u_mid,zeros(2*Nx-2,1));
            F_u = M'*Phi+F;
            du = Phi_mid-Phi;
            dF = beta*(F_u_mid-F_u);
            rho = norm(dF)/norm(du);
        end
        duF = du+dF;
        r1 = dot(du,du);
        r2 = dot(duF,duF);
        alpha = r1/r2;
        Phi = Phi_mid-alpha*rho0*duF;
        F_u = M'*Phi+F;
        tol = norm(Phi-max(Phi-F_u,zeros(2*Nx-2,1)),inf);
        if (rho<mu)
           beta = beta*rho0; 
        end
        k = k+1;
    end     
    Phi_end = Phi;
    
end
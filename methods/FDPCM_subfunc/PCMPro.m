function Phi_end = PCMPro(eps,nu,mu,vrho,M,F,Phi_start)
    % pcm�㷨���ĳ���

    % Nx�ǿռ䷽���ʷֶ���
    % eps����ֹ�������
    % nu��mu��vrho�Ǹ�����(0,1)֮���ϵ��
    % M�Ǹ�����(2Nx-2)*(2Nx-2)��ϵ������
    % F�Ǹ�����(2Nx-2)*1��������
    % Phi_start�Ǹò������ĳ�ʼֵ
    
    len = size(F,1);
    Phi = Phi_start;
    
    beta = 1;
    k = 0;
    F_u = M'*Phi+F;
    tol = norm(Phi-max(Phi-F_u,zeros(len,1)),inf);
    while (tol>eps&&k<10000)
        % ͶӰ
        Phi_mid = Phi;
        F_u_mid = F_u;
        Phi = max(Phi_mid-beta.*F_u_mid,0);
        % ���з���
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
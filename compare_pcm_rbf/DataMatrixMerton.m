function [Vtrue1,Vtrue2] = DataMatrixMerton(K,sigma1,sigma2,lam1,lam2,gamma,mun,rho0,eps,L,A,x0,t0,Nx,Nt,dx,dt,...
    mu,delt,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2)
    % 求解出分划网格各点V函数对应的值，存放进Vtrue矩阵

    
    % Phi_pre用于存放n层向量(已知)
    % Phi_cur用于存放n+1层向量(待求)
    % Vtrue存放还原后的数值
    Vtrue1 = zeros(Nt+1,Nx+1);
    Vtrue2 = zeros(Nt+1,Nx+1);
    
    x = linspace(-L,L,Nx+1);
    S = exp(x);
    
    % 迭代循环赋初值
    % 形成矩阵Phi_star
    x = linspace(-L,L,Nx+1);
    w1_star = exp(-nM1*x).*max(K-exp(x),0);
    Phi1_star = w1_star(2:1:Nx).';
    w2_star = exp(-nM2*x).*max(K-exp(x),0);
    Phi2_star = w2_star(2:1:Nx).';
    Phi_pre = [Phi1_star;Phi2_star];  
    
    % 赋值R矩阵
    [Merton1,Merton2] = AssembleRMerton(Nx,delt,mu_star1,mu_star2,dx,C1,C2);
    [B, d1t0, d2t0, F0t0, Phi_star_t0] = AssembleBasicMatMerton(dt,dx,Nx,sigma1,sigma2,nM1,nM2,mM1,mM2,A,lam1,lam2,L,K,...
    mu_star1,mu_star2,mu,delt);
    % 再给矩阵其余部分赋值
    for i = 2:Nt+1
        [F, Phi_star] = AssemblecurrMerton(i-1,dt,Nx,mM1,mM2,lam1,lam2,Phi_pre,Merton1,Merton2,B,d1t0,d2t0,F0t0,Phi_star_t0);
        Phi_cur = PCMPro(Nx,eps,gamma,mun,rho0,B,F,Phi_pre-Phi_star);
        Vtrue1(i,2:Nx) = Phi_cur(1:Nx-1).'+Phi_star(1:Nx-1).';
        Vtrue2(i,2:Nx) = Phi_cur(Nx:2*Nx-2).'+Phi_star(Nx:2*Nx-2).';
        [V1,V2] = RevValueMerton(i-2,Nx,Phi_cur,Phi_star,mM1,mM2,nM1,nM2,t0,x0,dt,dx);
        Vtrue1(i,2:Nx) = V1.';
        Vtrue2(i,2:Nx) = V2.';
        Phi_pre = Phi_cur+Phi_star;
    end
    
    % 还原矩阵边界
    % 上侧
    for i = 2:Nx
        Vtrue1(1,i) = V_star(K,S(i));
        Vtrue2(1,i) = V_star(K,S(i));
    end
    % 右侧
    for i = 1:Nt+1
        Vtrue1(i,Nx+1) = 0;
        Vtrue2(i,Nx+1) = 0;
    end
    % 左侧
    for i = 1:Nt+1
        Vtrue1(i,1) = V_star(K,exp(-L));
        Vtrue2(i,1) = V_star(K,exp(-L));
    end
    
end
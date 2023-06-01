function [Vtrue1,Vtrue2] = SolveRBF(Lmat,Umat,A1,L,Nx,Nt,K)
    % 求解出分划网格各点V函数对应的值，存放进Vtrue矩阵
    
    % Phi_pre用于存放n层向量(已知)
    % Phi_cur用于存放n+1层向量(待求)
    % Vtrue存放最终数值
    Vtrue1 = zeros(Nt+1,Nx+1);
    Vtrue2 = zeros(Nt+1,Nx+1);
    
    x = linspace(-L,L,Nx+1);
    S = exp(x);
    
    % 矩阵边界赋值
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
    
    % 再给矩阵其余部分赋值
    for i = 2:Nt+1
        Vpre = [Vtrue1(i-1,:) Vtrue2(i-1,:)]';
        mid = Lmat\Vpre;
        rho = Umat\mid;
        rho1 = rho(1:Nx+1,1);
        rho2 = rho(Nx+2:end,1);
        for j = 2:Nx
            Vtrue1(i,j) = max(max(K-exp(x(j)),0),A1(j,:)*rho1);
            Vtrue2(i,j) = max(max(K-exp(x(j)),0),A1(j,:)*rho2);
        end
    end
    
end
function [Vtrue1,Vtrue2] = SolveRBF(Lmat,Umat,A1,L,Nx,Nt,K)
    % Solve the values of the V function at each grid point and store them in the Vtrue matrix
    
    % Phi_pre is used to store the nth layer vector (known)
    % Phi_cur is used to store the n+1 th layer vector (to be solved)
    % Vtrue stores the final numerical values
    Vtrue1 = zeros(Nt+1, Nx+1);
    Vtrue2 = zeros(Nt+1, Nx+1);
    
    x = linspace(-L, L, Nx+1);
    S = exp(x);
    
    % Assign values to the matrix boundaries
    % Top side
    for i = 2:Nx
        Vtrue1(1, i) = V_star(K, S(i));
        Vtrue2(1, i) = V_star(K, S(i));
    end
    % Right side
    for i = 1:Nt+1
        Vtrue1(i, Nx+1) = 0;
        Vtrue2(i, Nx+1) = 0;
    end
    % Left side
    for i = 1:Nt+1
        Vtrue1(i, 1) = V_star(K, exp(-L));
        Vtrue2(i, 1) = V_star(K, exp(-L));
    end
    
    % Assign values to the remaining part of the matrix
    for i = 2:Nt+1
        Vpre = [Vtrue1(i-1, :) Vtrue2(i-1, :)]';
        mid = Lmat \ Vpre;
        rho = Umat \ mid;
        rho1 = rho(1:Nx+1, 1);
        rho2 = rho(Nx+2:end, 1);
        for j = 2:Nx
            Vtrue1(i, j) = max(max(K - exp(x(j)), 0), A1(j, :) * rho1);
            Vtrue2(i, j) = max(max(K - exp(x(j)), 0), A1(j, :) * rho2);
        end
    end
    
end
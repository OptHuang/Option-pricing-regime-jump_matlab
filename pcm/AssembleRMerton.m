function [Merton1,Merton2] = AssembleRMerton(Nx,delt,mu_star1,mu_star2,dx,C1,C2)
    
    Merton1 = zeros(Nx-1);
    Merton2 = zeros(Nx-1);
    
    for j = 1:Nx-1
        for k = 1:Nx-1
            Merton1(j,k) = f(delt,mu_star1,dx,j,k);
            Merton2(j,k) = f(delt,mu_star2,dx,j,k);
        end
    end
    
    Merton1 = C1*dx*Merton1;
    Merton2 = C2*dx*Merton2;
    
end
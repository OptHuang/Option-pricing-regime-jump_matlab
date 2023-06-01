function [F, Phi_star] = AssemblecurrMerton(n,dt,Nx,mM1,mM2,lam1,lam2,Phi_pre,Merton1,Merton2,B,d1t0,d2t0,F0t0,Phi_star_t0)

    %%%%%%%%%%%%%%形成向量F, F=F0+(A-R)*Phi(n)+B*Phi_star
    F = zeros(2*Nx-2,1);
    
    % Phi_star
    Phi_star = zeros(2*Nx-2,1);
    par1 = exp(-mM1*(n-1)*dt);
    par2 = exp(-mM2*(n-1)*dt);
    Phi_star(1:Nx-1,1) = par1*Phi_star_t0(1:Nx-1,1);
    Phi_star(Nx:end,1) = par2*Phi_star_t0(Nx:end,1);

    % A*Phi(n) - R*Phi(n) + B*Phi_star
    C1 = exp((mM2-mM1)*(n-1)*dt);
    F(1:Nx-1,1) = -Phi_pre(1:Nx-1,1) + C1*d1t0.*Phi_pre(Nx:end,1) - dt*lam1*Merton1*Phi_pre(1:Nx-1,1);
    F(Nx:end,1) = 1/C1*d2t0.*Phi_pre(1:Nx-1,1) - Phi_pre(Nx:end,1) - dt*lam2*Merton2*Phi_pre(Nx:end,1);
    F = F + B*Phi_star;
    
    % F0
    F(1:Nx-1,1) = F(1:Nx-1,1) + par1*F0t0(1:Nx-1,1);
    F(Nx:end,1) = F(Nx:end,1) + par2*F0t0(Nx:end,1);

end
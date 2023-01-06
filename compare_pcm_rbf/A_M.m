function res = A_M(K,L,mu,delt,nM1,nM2,mM1,mM2,dx,dt,i,j,n)
    
    switch i
        case 1
            res = exp(-nM1*(-L+j*dx)-mM1*n*dt)*(K*exp(-delt^2*nM1^2)*normcdf(-j*dx,mu,delt)-...
                exp(mu+0.5*delt^2-L+j*dx-delt^2*nM1^2)*normcdf(-j*dx,mu+delt^2,delt));
        case 2
            res = exp(-nM2*(-L+j*dx)-mM2*n*dt)*(K*exp(-delt^2*nM2^2)*normcdf(-j*dx,mu,delt)-...
                exp(mu+0.5*delt^2-L+j*dx-delt^2*nM2^2)*normcdf(-j*dx,mu+delt^2,delt));
    end
end
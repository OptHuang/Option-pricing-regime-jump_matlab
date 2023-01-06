function res = A_K(K,L,q,a2,nK1,nK2,mK1,mK2,dx,dt,i,j,n)

    switch i
        case 1
            res = q*exp(nK1*L-nK1*j*dx-a2*j*dx-mK1*n*dt)*(K-a2/(1+a2)*exp(-L));
        case 2
            res = q*exp(nK2*L-nK2*j*dx-a2*j*dx-mK2*n*dt)*(K-a2/(1+a2)*exp(-L));
    end
end
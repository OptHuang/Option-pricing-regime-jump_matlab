function val = f(delt,mu_star,dx,j,k)
    
    val = exp(-((k-j)*dx-mu_star)^2/(2*delt^2));
    
end
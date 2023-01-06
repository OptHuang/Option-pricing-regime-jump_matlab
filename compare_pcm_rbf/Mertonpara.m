function [mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2] = Mertonpara(sigma1,sigma2,r1,r2,d1,d2,lam1,lam2,A)
    % Merton型的参数
    
    mu = -0.025;
    delt = sqrt(0.05);  
    kappaM = exp(mu+delt^2/2)-1;
    nM1 = 1/2+1/(sigma1^2)*(-r1+d1+lam1*kappaM);
    nM2 = 1/2+1/(sigma2^2)*(-r2+d2+lam2*kappaM);
    mM1 = -(1/2*sigma1^2-r1+d1+lam1*kappaM)*nM1+1/2*sigma1^2*nM1^2-r1-lam1+A(1,1);
    mM2 = -(1/2*sigma2^2-r2+d2+lam2*kappaM)*nM2+1/2*sigma2^2*nM2^2-r2-lam2+A(2,2);   
    m1 = 1/2*sigma1^2*nM1^2-r1-lam1+A(1,1);
    mu_star1 = mu+delt^2*nM1;
    mu_star2 = mu+delt^2*nM2;
    C1 = exp(mu*nM1-1/2*delt^2*nM1^2)/delt/sqrt(2*pi);
    C2 = exp(mu*nM2-1/2*delt^2*nM2^2)/delt/sqrt(2*pi);

end
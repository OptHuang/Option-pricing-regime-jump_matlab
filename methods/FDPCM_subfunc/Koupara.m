function [p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = Koupara(sigma1,sigma2,r1,r2,d1,d2,lam1,lam2,A)
    % 输入Kou型参数
    
    p = 0.3445;  
    q = 1-p;
    a1 = 3.0465;
    a2 = 3.0775;
    kappaK = p*a1/(a1-1)+(1-p)*a2/(a2+1)-1;
    nK1 = 1/2+1/(sigma1^2)*(-r1+d1+lam1*kappaK);
    nK2 = 1/2+1/(sigma2^2)*(-r2+d2+lam2*kappaK);
    mK1 = -(1/2*sigma1^2-r1+d1+lam1*kappaK)*nK1+1/2*sigma1^2*nK1^2-r1-lam1+A(1,1);
    mK2 = -(1/2*sigma2^2-r2+d2+lam2*kappaK)*nK2+1/2*sigma2^2*nK2^2-r2-lam2+A(2,2);  
    
end
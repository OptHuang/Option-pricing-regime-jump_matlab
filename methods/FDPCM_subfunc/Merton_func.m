function [y,deriv_y] = Merton_func(sigma,r,d,lam,kappaM,mu,delt,x)
    % Merton型的特征函数, 返回特征函数的值和导数值
    
    y = 0.5*sigma^2*x^2+(r-d-lam*kappaM-0.5*sigma^2)*x-(r+lam)+lam*exp(mu*x+0.5*delt^2*x^2);
    deriv_y = sigma^2*x+(r-d-lam*kappaM-0.5*sigma^2)+lam*(mu+sigma^2*x)*exp(mu*x+0.5*delt^2*x^2);

end
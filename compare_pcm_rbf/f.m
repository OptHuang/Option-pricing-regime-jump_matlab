function val = f(delt,mu_star,dx,j,k)
    % 计算Merton类型中的fij(k)函数
    
    val = exp(-((k-j)*dx-mu_star)^2/(2*delt^2));
    
end
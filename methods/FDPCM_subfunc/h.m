function val = h(qp,j,k,nK,dx,p,q,a1,a2)
    % 计算Kou型中的Qi*hij^Q(k)或Pi*hij^P(k)的值
    % qp:0表示Q,1表示P
    % n对应region输入n1或n2
    
    switch qp
        case 0
            par = q*a2;
            val = par*exp((nK+a2)*(k-j)*dx);
        case 1
            par = p*a1;
            val = par*exp((nK-a1)*(k-j)*dx);
    end
    
end
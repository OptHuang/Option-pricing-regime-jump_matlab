function val = h(qp,j,k,nK,dx,p,q,a1,a2)
    % ����Kou���е�Qi*hij^Q(k)��Pi*hij^P(k)��ֵ
    % qp:0��ʾQ,1��ʾP
    % n��Ӧregion����n1��n2
    
    switch qp
        case 0
            par = q*a2;
            val = par*exp((nK+a2)*(k-j)*dx);
        case 1
            par = p*a1;
            val = par*exp((nK-a1)*(k-j)*dx);
    end
    
end
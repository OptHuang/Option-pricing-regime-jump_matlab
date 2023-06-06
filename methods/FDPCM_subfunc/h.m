function val = h(qp,j,k,nK,dx,p,q,a1,a2)
    % Calculate the value of Qihij^Q(k) or Pihij^P(k) in Kou type model
    % qp: 0 represents Q, 1 represents P
    % n corresponds to region input n1 or n2
    
    switch qp
        case 0
            par = q*a2;
            val = par*exp((nK+a2)*(k-j)*dx);
        case 1
            par = p*a1;
            val = par*exp((nK-a1)*(k-j)*dx);
    end
    
end
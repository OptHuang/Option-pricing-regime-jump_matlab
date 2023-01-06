function [Kou1,Kou2] = AssembleRKou(Nx,nK1,nK2,dx,p,q,a1,a2)

    Kou1 = zeros(Nx-1);
    Kou2 = zeros(Nx-1);
    
    for j = 1:Nx-1
        Kou1(j,j) = 1/2*(h(0,j,j-1,nK1,dx,p,q,a1,a2)+h(1,j,j,nK1,dx,p,q,a1,a2));
        Kou2(j,j) = 1/2*(h(0,j,j-1,nK2,dx,p,q,a1,a2)+h(1,j,j,nK2,dx,p,q,a1,a2));
    end
    
    for j = 2:Nx-1
        for k = 1:j-1
            Kou1(j,k) = h(0,j,k,nK1,dx,p,q,a1,a2);
            Kou2(j,k) = h(0,j,k,nK2,dx,p,q,a1,a2);
        end
    end
    
    for j = 1:Nx-2
        for k = j+1:Nx-1
            Kou1(j,k) = h(1,j,k,nK1,dx,p,q,a1,a2);
            Kou2(j,k) = h(1,j,k,nK2,dx,p,q,a1,a2);
        end
    end
    
    Kou1 = dx*Kou1;
    Kou2 = dx*Kou2;
    
end
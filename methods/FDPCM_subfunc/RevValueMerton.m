function [V1,V2] = RevValueMerton(n,Nx,Phi_cur,Phi_star,mM1,mM2,nM1,nM2,t0,x0,dt,dx)
    % This subfunction is used to restore the solution of the n+1 layer obtained by the pcm algorithm to the values of the corresponding points of the required function V
    
    % n refers to the nth layer in the time partition, Nx is the number of space partition segments
    % Phi_cur is the (2Nx-2)*1 column vector obtained directly by pcm
    % G is the corresponding (2Nx-2)*1 column vector for this layer
    % V1, V2 are the values of the Value function under regimes 1 and 2 at the selected points
    % mM1, mM2, nM1, nM2 are the parameters in the transformation calculated by AssembleMatrix
    % t0, x0, dt, dx are the starting point of time, the starting point of space, the length of each segment of the time partition, and the length of each segment of the space partition
    
    Phi = Phi_cur+Phi_star;
    W1 = Phi(1:Nx-1);
    W2 = Phi(Nx:2*Nx-2);
    
    c1 = zeros(1,Nx-1);
    c2 = zeros(1,Nx-1);
    for i =1:Nx-1
        c1(i) = exp(mM1*(t0+(n+1)*dt)+nM1*(x0+i*dx));
        c2(i) = exp(mM2*(t0+(n+1)*dt)+nM2*(x0+i*dx));
    end
    reC1 = diag(c1);
    reC2 = diag(c2);
    
    V1 = reC1*W1;
    V2 = reC2*W2;

end
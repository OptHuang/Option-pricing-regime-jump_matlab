function [V1,V2] = RevValueKou(n,Nx,Phi_cur,Phi_star,mK1,mK2,nK1,nK2,t0,x0,dt,dx)
    % 这个子函数用于将pcm算法求得的第n+1层的解还原成所需要的函数V所对应各点的值
    
    % n是指在时间分划第n层，Nx是空间分划段数
    % Phi_cur是pcm直接求得的(2Nx-2)*1的列向量
    % G是该层对应的(2Nx-2)*1的列向量
    % V1,V2就是制度1,2下分别对应的Value函数在所取的点处的值
    % m1,m2,n1,n2是AssembleMatrix计算出来的转化中的参数
    % t0,x0,dt,dx分别是时间起点，空间起点，时间分划每段长度，空间分划每段长度
    
    Phi = Phi_cur+Phi_star;
    W1 = Phi(1:Nx-1);
    W2 = Phi(Nx:2*Nx-2);
    
    c1 = zeros(1,Nx-1);
    c2 = zeros(1,Nx-1);
    for i =1:Nx-1
        c1(i) = exp(mK1*(t0+(n+1)*dt)+nK1*(x0+i*dx));
        c2(i) = exp(mK2*(t0+(n+1)*dt)+nK2*(x0+i*dx));
    end
    reC1 = diag(c1);
    reC2 = diag(c2);
    
    V1 = reC1*W1;
    V2 = reC2*W2;

end
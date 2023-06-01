function [V1,V2] = RevValueMerton(n,Nx,Phi_cur,Phi_star,mM1,mM2,nM1,nM2,t0,x0,dt,dx)
    % ����Ӻ������ڽ�pcm�㷨��õĵ�n+1��Ľ⻹ԭ������Ҫ�ĺ���V����Ӧ�����ֵ
    
    % n��ָ��ʱ��ֻ���n�㣬Nx�ǿռ�ֻ�����
    % Phi_cur��pcmֱ����õ�(2Nx-2)*1��������
    % G�Ǹò��Ӧ��(2Nx-2)*1��������
    % V1,V2�����ƶ�1,2�·ֱ��Ӧ��Value��������ȡ�ĵ㴦��ֵ
    % mM1,mM2,nM1,nM2��AssembleMatrix���������ת���еĲ���
    % t0,x0,dt,dx�ֱ���ʱ����㣬�ռ���㣬ʱ��ֻ�ÿ�γ��ȣ��ռ�ֻ�ÿ�γ���
    
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
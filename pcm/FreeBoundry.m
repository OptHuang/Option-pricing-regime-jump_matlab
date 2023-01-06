function [fb1_num,fb2_num,fb1,fb2] = FreeBoundry(Utrue1,Utrue2,K,L,Nx,Nt)
    % 返回的是一个向量
    % 每一个分量表示自由边界在对应时间层在Utrue中的位置
    
    fb1_num = zeros(Nt+1,1);
    fb2_num = zeros(Nt+1,1);
    fb1 = zeros(Nt,1);
    fb2 = zeros(Nt,1);
    x = linspace(-L,L,Nx+1);
    s = exp(x);
    
    eps = 1e-8;
    
    for i = 2:Nt+1
        fb1_num(i) = min(find(abs(Utrue1(i,:)-max(K-s,0))>eps));
        fb1(i) = Utrue1(i,fb1_num(i));
    end

    
    for i = 2:Nt+1
        fb2_num(i) = min(find(abs(Utrue2(i,:)-max(K-s,0))>eps));
        fb2(i) = Utrue2(i,fb2_num(i));
    end
     
end
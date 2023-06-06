function [fb1_num,fb2_num,fb1,fb2] = FreeBoundry(Utrue1,Utrue2,K,L,Nx,Nt)
    
    fb1_num = zeros(Nt+1,1);
    fb2_num = zeros(Nt+1,1);
    fb1 = zeros(Nt,1);
    fb2 = zeros(Nt,1);
    x = linspace(-L,L,Nx+1);
    s = exp(x);
    
    eps = 5e-5;
    
    for i = 2:Nt+1
        temp_result = find(abs(Utrue1(i,:)-max(K-s,0))>eps, 1);
        if isempty(temp_result)
            fb1_num(i) = find(max(K-s,0)==0, 1);
        else
            fb1_num(i) = temp_result;
        end
        fb1(i) = Utrue1(i,fb1_num(i));
    end
    
    for i = 2:Nt+1
        temp_result = find(abs(Utrue2(i,:)-max(K-s,0))>eps, 1);
        if isempty(temp_result)
            fb2_num(i) = find(max(K-s,0)==0, 1);
        else
            fb2_num(i) = temp_result;
        end
        fb2(i) = Utrue2(i,fb2_num(i));
    end
     
end
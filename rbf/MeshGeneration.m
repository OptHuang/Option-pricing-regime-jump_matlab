function [S,M,t,s] = MeshGeneration(T,L,Nx,Nt)
    % 二维网格剖分函数(为绘图做铺垫)
    
    % T,L：剖分矩形区域的边界
    % Nx：空间方向的剖分结点段数
    % Nt：时间方向的剖分结点段数
    
    % 找到原始坐标轴上与差分法所取的点相对应的点
    t = linspace(0,T,Nt+1);
    x = linspace(-L,L,Nx+1); 
    s = exp(x);
    
    [S,M] = meshgrid(s,t);
    
end
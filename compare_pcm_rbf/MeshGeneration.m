function [S,M,t,s] = MeshGeneration(T,L,Nx,Nt)
    % ��ά�����ʷֺ���(Ϊ��ͼ���̵�)
    
    % T,L���ʷ־�������ı߽�
    % Nx���ռ䷽����ʷֽ�����
    % Nt��ʱ�䷽����ʷֽ�����
    
    % �ҵ�ԭʼ�����������ַ���ȡ�ĵ����Ӧ�ĵ�
    t = linspace(0,T,Nt+1);
    x = linspace(-L,L,Nx+1);
    s = zeros(1,Nx+2);
    s(1) = 0;
    s(2:end) = exp(x);
    
    [S,M] = meshgrid(s,t);
    
end
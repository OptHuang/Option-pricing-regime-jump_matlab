function Plot3D(method,Nx,Nt,problem_paras,options)
    
    currentFolder = fileparts(mfilename('fullpath'));
    testSubFolder = fullfile(currentFolder, '..', 'methods');
    addpath(testSubFolder);

    switch method
        case {"FDPCM"}
            if nargin < 5
                [~,~,full_matrix1,full_matrix2,fb1,fb2,L] = FDPCM(Nx,Nt,problem_paras);
            else
                [~,~,full_matrix1,full_matrix2,fb1,fb2,L] = FDPCM(Nx,Nt,problem_paras,options);
            end
        case {"RBCM"}
            if nargin < 5
                [~,~,full_matrix1,full_matrix2,fb1,fb2,L] = RBCM(Nx,Nt,problem_paras);
            else
                [~,~,full_matrix1,full_matrix2,fb1,fb2,L] = RBCM(Nx,Nt,problem_paras,options);
            end
    end
    
    t = linspace(0,problem_paras.T,Nt+1);
    x = linspace(-L,L,Nx+1); 
    s = exp(x);
    [S,M] = meshgrid(s,t);

    full_matrix1 = flipud(full_matrix1);
    full_matrix2 = flipud(full_matrix2);
    fb1 = flipud(fb1);
    fb2 = flipud(fb2);
    
    subplot(1,2,1)
    mesh(S,M,full_matrix1)
    hold on
    plot3(problem_paras.K-fb1,t,fb1,'k-','LineWidth',3)
    hold on
    plot3(problem_paras.K-fb1,t,zeros(size(t)),'k-','LineWidth',3)
    axis([0 1.2*exp(L) 0 problem_paras.T 0 1.2*problem_paras.K])
    xlabel('S')
    ylabel('T')
    zlabel('P')
    switch problem_paras.type
        case {"Merton"}
            title('Example1: option value under regime 1')
        case {"Kou"}
            title('Example2: option value under regime 1')
    end
    
    subplot(1,2,2)
    mesh(S,M,full_matrix2)
    hold on
    plot3(problem_paras.K-fb2,t,fb2,'k-','LineWidth',3)
    hold on
    plot3(problem_paras.K-fb1,t,zeros(size(t)),'k-','LineWidth',3)
    axis([0 1.2*exp(L) 0 problem_paras.T 0 1.2*problem_paras.K])
    xlabel('S')
    ylabel('T')
    zlabel('P')
    switch problem_paras.type
        case {"Merton"}
            title('Example1: option value under regime 2')
        case {"Kou"}
            title('Example2: option value under regime 2')
    end
    
            
end
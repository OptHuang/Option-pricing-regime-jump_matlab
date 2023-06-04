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

    full_matrix1 = [problem_paras.K*ones(Nt+1,1), full_matrix1];
    full_matrix2 = [problem_paras.K*ones(Nt+1,1), full_matrix2];
    
    t = linspace(0,problem_paras.T,Nt+1);
    x = linspace(-L,L,Nx+1); 
    s = exp(x);
    s = [0, s];
    [S,M] = meshgrid(s,t);

    full_matrix1 = flipud(full_matrix1);
    full_matrix2 = flipud(full_matrix2);
    fb1 = flipud(fb1);
    fb2 = flipud(fb2);
    
    subplot(1,2,1)
    mesh(S,M,full_matrix1)
    hold on
    plot3(problem_paras.K-fb1,t,fb1,'k-','LineWidth',3)
%     hold on
%     plot3(problem_paras.K-fb1,t,zeros(size(t)),'k-','LineWidth',3)
    axis([0 exp(L) 0 problem_paras.T 0 problem_paras.K])
    xlabel('$S$','Interpreter','latex','FontSize',10)
    ylabel('$T$','Interpreter','latex','FontSize',10)
    zlabel('$P$','Interpreter','latex','FontSize',10)
    switch problem_paras.type
        case {"Merton"}
            title('Option price under regime 1 in Example 1','Interpreter','latex','FontSize',12)
        case {"Kou"}
            title('Option price under regime 1 in Example 2','Interpreter','latex','FontSize',12)
    end
    
    subplot(1,2,2)
    mesh(S,M,full_matrix2)
    hold on
    plot3(problem_paras.K-fb2,t,fb2,'k-','LineWidth',3)
%     hold on
%     plot3(problem_paras.K-fb1,t,zeros(size(t)),'k-','LineWidth',3)
    axis([0 exp(L) 0 problem_paras.T 0 problem_paras.K])
    xlabel('$S$','Interpreter','latex','FontSize',10)
    ylabel('$T$','Interpreter','latex','FontSize',10)
    zlabel('$P$','Interpreter','latex','FontSize',10)
    switch problem_paras.type
        case {"Merton"}
            title('Option price under regime 2 in Example 1','Interpreter','latex','FontSize',12)
        case {"Kou"}
            title('Option price under regime 2 in Example 2','Interpreter','latex','FontSize',12)
    end
              
end
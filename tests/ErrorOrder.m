function ErrorOrder(fixwhich,npoint,problem_paras,Nx_true,mesh_ratio_true,options)

    % fixwhich: the input parameter to decide which variable to be
    % fixed; if fixwhich = "spatial", then fix variable "Nx"; if fixwhich =
    % "temporal", then fix variable "Nt"; if fixwhich = "ratio", then fix
    % variable mesh_ratio.
    %
    % npoint: the number of points to draw in the picture.

    currentFolder = fileparts(mfilename('fullpath'));
    testSubFolder = fullfile(currentFolder, '..', 'methods');
    addpath(testSubFolder);
    TruncationTechFolder = fullfile(currentFolder, '..', 'methods/FDPCM_subfunc');
    addpath(TruncationTechFolder);

    folder_path = '../data_TrueSol';
    file_pattern = sprintf('type=%s_ratio=%s_Nx=%d_Nt=*_d=%s_lam1=%s_lam2=%s.mat',...
        problem_paras.type, strrep(sprintf('%.2f', mesh_ratio_true), '.', '-'), Nx_true, ...
        strrep(sprintf('%.3f', problem_paras.d1), '.', '-'), ...
        strrep(sprintf('%.2f', problem_paras.lam1), '.', '-'), ...
        strrep(sprintf('%.2f', problem_paras.lam2), '.', '-'));
    
    files = dir(fullfile(folder_path, file_pattern));
    
    if isempty(files)
        fprintf('No matching files found.\n');
    else
        truedata = fullfile(folder_path, files(1).name);
        fprintf('Loading file: %s\n', truedata);
        load(truedata);
    end

    % After loading the truedata file
    pattern = 'type=.*_ratio=.*_Nx=.*_Nt=(?<Nt>\d+)_d=.*_lam1=.*_lam2=.*\.mat';
    tokens = regexp(files(1).name, pattern, 'names');
    Nt_true = str2double(tokens.Nt);



%     load("Vkoutrue10000512.mat")
%     Nt_true = 10000;
%     TrueSol_last_row1 = Vkoutrue1(end,:);
%     TrueSol_last_row2 = Vkoutrue2(end,:);



    Error1 = zeros(1,npoint);
    Error2 = zeros(1,npoint);
    N0 = zeros(1,npoint);

    if nargin < 6
        switch fixwhich
            case {"spatial"}
                Nt0 = 500;
                Nx = Nx_true;
    
                imax = floor(log2(Nt_true/Nt0));
                if npoint > imax
                    fprintf("npoint is too big!")
                    return;
                end
    
                for i = 1:npoint
                    Nt = Nt0*2^(i-1);
                    [last_row_pcm1,last_row_pcm2,~,~,~,~,L] = FDPCM(Nx,Nt,problem_paras);
                    intF1 = last_row_pcm1 - TrueSol_last_row1;
                    intF2 = last_row_pcm2 - TrueSol_last_row2;
                    Error1(i) = sqrt(intF1*intF1'*2*L/Nx);
                    Error2(i) = sqrt(intF2*intF2'*2*L/Nx);
                end
    
                for i = 1:npoint
                    N0(i) = Nt0*2^(i-1);
                end
            
                k = 1;
                Standard = -k*log(N0)+1;
    
                max_xlabel = max(log(N0));
                min_xlabel = min(log(N0));
                len_xlabel = max_xlabel - min_xlabel;
                left_xlabel = -0.1*len_xlabel + min_xlabel;
                right_xlabel = 0.1*len_xlabel + max_xlabel;
                max_ylabel = max([max(Standard), max(log(Error1)), max(log(Error2))]);
                min_ylabel = min([min(Standard), min(log(Error1)), min(log(Error2))]);
                len_ylabel = max_ylabel - min_ylabel;
                lower_ylabel = -0.1*len_ylabel + min_ylabel;
                upper_ylabel = 0.1*len_ylabel + max_ylabel;

                plot(log(N0),Standard,'-b',log(N0),log(Error1),'-ro',log(N0),log(Error2),'-g*')
                xlabel('$\log(N_t)$','Interpreter','latex','FontSize',16)
                ylabel('log(error)','Interpreter','latex','FontSize',16)
                axis([left_xlabel, right_xlabel, lower_ylabel, upper_ylabel])
                legend('order 1','log(error)--regime 1','log(error)--regime 2','Interpreter','latex','FontSize',12)
                switch problem_paras.type
                    case {"Merton"}
                        title('Convergence order in Example 1','Interpreter','latex','FontSize',22)
                    case {"Kou"}
                        title('Convergence order in Example 2','Interpreter','latex','FontSize',22)
                end
            case {"temporal"}
                Nx0 = 16;
                Nt = Nt_true;
        
                imax = log2(Nx_true/Nx0);
                if npoint > imax
                    fprintf("npoint is too big!")
                    return;
                end
    
                for i = 1:npoint
                    Nx = Nx0*2^(i-1);
                    [last_row_pcm1,last_row_pcm2,~,~,~,~,L] = FDPCM(Nx,Nt,problem_paras);
                    mult = Nx_true/Nx0/2^(i-1);
                    True_piece1 = zeros(size(last_row_pcm1));
                    True_piece2 = zeros(size(last_row_pcm1));
                    for j = 1:Nx0*2^(i-1)
                        True_piece1(j) = TrueSol_last_row1(1+(j-1)*mult);
                        True_piece2(j) = TrueSol_last_row2(1+(j-1)*mult);
                    end
                    intF1 = last_row_pcm1 - True_piece1;
                    intF2 = last_row_pcm2 - True_piece2;
                    Error1(i) = sqrt(intF1*intF1'*2*L/(Nx0*2^(i-1)));
                    Error2(i) = sqrt(intF2*intF2'*2*L/(Nx0*2^(i-1)));
                end
                
                for i = 1:npoint
                    N0(i) = Nx0*2^(i-1);
                end
            
                k = 2;
                Standard = -k*log(N0)+1;

                max_xlabel = max(log(N0));
                min_xlabel = min(log(N0));
                len_xlabel = max_xlabel - min_xlabel;
                left_xlabel = -0.1*len_xlabel + min_xlabel;
                right_xlabel = 0.1*len_xlabel + max_xlabel;
                max_ylabel = max([max(Standard), max(log(Error1)), max(log(Error2))]);
                min_ylabel = min([min(Standard), min(log(Error1)), min(log(Error2))]);
                len_ylabel = max_ylabel - min_ylabel;
                lower_ylabel = -0.1*len_ylabel + min_ylabel;
                upper_ylabel = 0.1*len_ylabel + max_ylabel;
    
                plot(log(N0),Standard,'-b',log(N0),log(Error1),'-ro',log(N0),log(Error2),'-g*')
                xlabel('$\log(N_x)$','Interpreter','latex','FontSize',16)
                ylabel('log(error)','Interpreter','latex','FontSize',16)
                axis([left_xlabel, right_xlabel, lower_ylabel, upper_ylabel])
                legend('order 2','log(error)--regime 1','log(error)--regime 2','Interpreter','latex','FontSize',12)
                switch problem_paras.type
                    case {"Merton"}
                        title('Convergence order in Example 1','Interpreter','latex','FontSize',22)
                    case {"Kou"}
                        title('Convergence order in Example 2','Interpreter','latex','FontSize',22)
                end
            case {"ratio"}
                Nx0 = 16;
                mesh_ratio = mesh_ratio_true;
        
                imax = log2(Nx_true/Nx0);
                if npoint > imax
                    fprintf("npoint is too big!")
                    return;
                end
    
                switch problem_paras.type
                        case {"Merton"}
                            L = TruncationTech(problem_paras.K,0,problem_paras.sig1,...
                                problem_paras.sig2,problem_paras.r1,problem_paras.r2,...
                                problem_paras.d1,problem_paras.d2,problem_paras.lam1,...
                                problem_paras.lam2,problem_paras.A);
                        case {"Kou"}
                            L = TruncationTech(problem_paras.K,1,problem_paras.sig1,...
                                problem_paras.sig2,problem_paras.r1,problem_paras.r2,...
                                problem_paras.d1,problem_paras.d2,problem_paras.lam1,...
                                problem_paras.lam2,problem_paras.A);
                end
    
                for i = 1:npoint
                    Nx = Nx0*2^(i-1);
                    Nt = floor(problem_paras.T / ((2*L/Nx)^2/mesh_ratio));
    
                    [last_row_pcm1,last_row_pcm2,~,~,~,~,~] = FDPCM(Nx,Nt,problem_paras);
                    mult = Nx_true/Nx0/2^(i-1);
                    True_piece1 = zeros(size(last_row_pcm1));
                    True_piece2 = zeros(size(last_row_pcm1));
                    for j = 1:Nx0*2^(i-1)
                        True_piece1(j) = TrueSol_last_row1(1+(j-1)*mult);
                        True_piece2(j) = TrueSol_last_row2(1+(j-1)*mult);
                    end
                    intF1 = last_row_pcm1 - True_piece1;
                    intF2 = last_row_pcm2 - True_piece2;
                    Error1(i) = sqrt(intF1*intF1'*2*L/(Nx0*2^(i-1)));
                    Error2(i) = sqrt(intF2*intF2'*2*L/(Nx0*2^(i-1)));
                end
                
                for i = 1:npoint
                    N0(i) = Nx0*2^(i-1);
                end
            
                k = 2;
                Standard = -k*log(N0)+1;

                max_xlabel = max(log(N0));
                min_xlabel = min(log(N0));
                len_xlabel = max_xlabel - min_xlabel;
                left_xlabel = -0.1*len_xlabel + min_xlabel;
                right_xlabel = 0.1*len_xlabel + max_xlabel;
                max_ylabel = max([max(Standard), max(log(Error1)), max(log(Error2))]);
                min_ylabel = min([min(Standard), min(log(Error1)), min(log(Error2))]);
                len_ylabel = max_ylabel - min_ylabel;
                lower_ylabel = -0.1*len_ylabel + min_ylabel;
                upper_ylabel = 0.1*len_ylabel + max_ylabel;
    
                plot(log(N0),Standard,'-b',log(N0),log(Error1),'-ro',log(N0),log(Error2),'-g*')
                xlabel('$\log(N_x)$','Interpreter','latex','FontSize',16)
                ylabel('log(error)','Interpreter','latex','FontSize',16)
                axis([left_xlabel, right_xlabel, lower_ylabel, upper_ylabel])
                legend('order 2','log(error)--regime 1','log(error)--regime 2','Interpreter','latex','FontSize',12)
                switch problem_paras.type
                    case {"Merton"}
                        title('Convergence order in Example 1','Interpreter','latex','FontSize',22)
                    case {"Kou"}
                        title('Convergence order in Example 2','Interpreter','latex','FontSize',22)
                end
        end
    else
        switch fixwhich
            case {"spatial"}
                Nt0 = 500;
                Nx = Nx_true;
    
                imax = floor(log2(Nt_true/Nt0));
                if npoint > imax
                    fprintf("npoint is too big!")
                    return;
                end
    
                for i = 1:npoint
                    Nt = Nt0*2^(i-1);
                    [last_row_pcm1,last_row_pcm2,~,~,~,~,L] = FDPCM(Nx,Nt,problem_paras,options);
                    intF1 = last_row_pcm1 - TrueSol_last_row1;
                    intF2 = last_row_pcm2 - TrueSol_last_row2;
                    Error1(i) = sqrt(intF1*intF1'*2*L/Nx);
                    Error2(i) = sqrt(intF2*intF2'*2*L/Nx);
                end
    
                for i = 1:npoint
                    N0(i) = Nt0*2^(i-1);
                end
            
                k = 1;
                Standard = -k*log(N0)+1;

                max_xlabel = max(log(N0));
                min_xlabel = min(log(N0));
                len_xlabel = max_xlabel - min_xlabel;
                left_xlabel = -0.1*len_xlabel + min_xlabel;
                right_xlabel = 0.1*len_xlabel + max_xlabel;
                max_ylabel = max([max(Standard), max(log(Error1)), max(log(Error2))]);
                min_ylabel = min([min(Standard), min(log(Error1)), min(log(Error2))]);
                len_ylabel = max_ylabel - min_ylabel;
                lower_ylabel = -0.1*len_ylabel + min_ylabel;
                upper_ylabel = 0.1*len_ylabel + max_ylabel;
    
                plot(log(N0),Standard,'-b',log(N0),log(Error1),'-ro',log(N0),log(Error2),'-g*')
                xlabel('$\log(N_t)$','Interpreter','latex','FontSize',16)
                ylabel('log(error)','Interpreter','latex','FontSize',16)
                axis([left_xlabel, right_xlabel, lower_ylabel, upper_ylabel])
                legend('order 1','log(error)--regime 1','log(error)--regime 2','Interpreter','latex','FontSize',12)
                switch problem_paras.type
                    case {"Merton"}
                        title('Convergence order in Example 1','Interpreter','latex','FontSize',22)
                    case {"Kou"}
                        title('Convergence order in Example 2','Interpreter','latex','FontSize',22)
                end
            case {"temporal"}
                Nx0 = 16;
                Nt = Nt_true;
        
                imax = log2(Nx_true/Nx0);
                if npoint > imax
                    fprintf("npoint is too big!")
                    return;
                end
    
                for i = 1:npoint
                    Nx = Nx0*2^(i-1);
                    [last_row_pcm1,last_row_pcm2,~,~,~,~,L] = FDPCM(Nx,Nt,problem_paras,options);
                    mult = Nx_true/Nx0/2^(i-1);
                    True_piece1 = zeros(size(last_row_pcm1));
                    True_piece2 = zeros(size(last_row_pcm1));
                    for j = 1:Nx0*2^(i-1)
                        True_piece1(j) = TrueSol_last_row1(1+(j-1)*mult);
                        True_piece2(j) = TrueSol_last_row2(1+(j-1)*mult);
                    end
                    intF1 = last_row_pcm1 - True_piece1;
                    intF2 = last_row_pcm2 - True_piece2;
                    Error1(i) = sqrt(intF1*intF1'*2*L/(Nx0*2^(i-1)));
                    Error2(i) = sqrt(intF2*intF2'*2*L/(Nx0*2^(i-1)));
                end
                
                for i = 1:npoint
                    N0(i) = Nx0*2^(i-1);
                end
            
                k = 2;
                Standard = -k*log(N0)+1;

                max_xlabel = max(log(N0));
                min_xlabel = min(log(N0));
                len_xlabel = max_xlabel - min_xlabel;
                left_xlabel = -0.1*len_xlabel + min_xlabel;
                right_xlabel = 0.1*len_xlabel + max_xlabel;
                max_ylabel = max([max(Standard), max(log(Error1)), max(log(Error2))]);
                min_ylabel = min([min(Standard), min(log(Error1)), min(log(Error2))]);
                len_ylabel = max_ylabel - min_ylabel;
                lower_ylabel = -0.1*len_ylabel + min_ylabel;
                upper_ylabel = 0.1*len_ylabel + max_ylabel;
    
                plot(log(N0),Standard,'-b',log(N0),log(Error1),'-ro',log(N0),log(Error2),'-g*')
                xlabel('$\log(N_x)$','Interpreter','latex','FontSize',16)
                ylabel('log(error)','Interpreter','latex','FontSize',16)
                axis([left_xlabel, right_xlabel, lower_ylabel, upper_ylabel])
                legend('order 2','log(error)--regime 1','log(error)--regime 2','Interpreter','latex','FontSize',12)
                switch problem_paras.type
                    case {"Merton"}
                        title('Convergence order in Example 1','Interpreter','latex','FontSize',22)
                    case {"Kou"}
                        title('Convergence order in Example 2','Interpreter','latex','FontSize',22)
                end
            case {"ratio"}
                Nx0 = 16;
                mesh_ratio = mesh_ratio_true;
        
                imax = log2(Nx_true/Nx0);
                if npoint > imax
                    fprintf("npoint is too big!")
                    return;
                end
    
                switch problem_paras.type
                        case {"Merton"}
                            L = TruncationTech(problem_paras.K,0,problem_paras.sig1,...
                                problem_paras.sig2,problem_paras.r1,problem_paras.r2,...
                                problem_paras.d1,problem_paras.d2,problem_paras.lam1,...
                                problem_paras.lam2,problem_paras.A);
                        case {"Kou"}
                            L = TruncationTech(problem_paras.K,1,problem_paras.sig1,...
                                problem_paras.sig2,problem_paras.r1,problem_paras.r2,...
                                problem_paras.d1,problem_paras.d2,problem_paras.lam1,...
                                problem_paras.lam2,problem_paras.A);
                end
    
                for i = 1:npoint
                    Nx = Nx0*2^(i-1);
                    Nt = floor(problem_paras.T / ((2*L/Nx)^2/mesh_ratio));
    
                    [last_row_pcm1,last_row_pcm2,~,~,~,~,~] = FDPCM(Nx,Nt,problem_paras,options);
                    mult = Nx_true/Nx0/2^(i-1);
                    True_piece1 = zeros(size(last_row_pcm1));
                    True_piece2 = zeros(size(last_row_pcm1));
                    for j = 1:Nx0*2^(i-1)
                        True_piece1(j) = TrueSol_last_row1(1+(j-1)*mult);
                        True_piece2(j) = TrueSol_last_row2(1+(j-1)*mult);
                    end
                    intF1 = last_row_pcm1 - True_piece1;
                    intF2 = last_row_pcm2 - True_piece2;
                    Error1(i) = sqrt(intF1*intF1'*2*L/(Nx0*2^(i-1)));
                    Error2(i) = sqrt(intF2*intF2'*2*L/(Nx0*2^(i-1)));
                end
                
                for i = 1:npoint
                    N0(i) = Nx0*2^(i-1);
                end
            
                k = 2;
                Standard = -k*log(N0)+1;

                max_xlabel = max(log(N0));
                min_xlabel = min(log(N0));
                len_xlabel = max_xlabel - min_xlabel;
                left_xlabel = -0.1*len_xlabel + min_xlabel;
                right_xlabel = 0.1*len_xlabel + max_xlabel;
                max_ylabel = max([max(Standard), max(log(Error1)), max(log(Error2))]);
                min_ylabel = min([min(Standard), min(log(Error1)), min(log(Error2))]);
                len_ylabel = max_ylabel - min_ylabel;
                lower_ylabel = -0.1*len_ylabel + min_ylabel;
                upper_ylabel = 0.1*len_ylabel + max_ylabel;
    
                plot(log(N0),Standard,'-b',log(N0),log(Error1),'-ro',log(N0),log(Error2),'-g*')
                xlabel('$\log(N_x)$','Interpreter','latex','FontSize',16)
                ylabel('log(error)','Interpreter','latex','FontSize',16)
                axis([left_xlabel, right_xlabel, lower_ylabel, upper_ylabel])
                legend('order 2','log(error)--regime 1','log(error)--regime 2','Interpreter','latex','FontSize',12)
                switch problem_paras.type
                    case {"Merton"}
                        title('Convergence order in Example 1','Interpreter','latex','FontSize',22)
                    case {"Kou"}
                        title('Convergence order in Example 2','Interpreter','latex','FontSize',22)
                end
        end
    end
    





end
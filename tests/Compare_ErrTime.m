function Compare_ErrTime(Nx,Nt,problem_paras,Nx_true,mesh_ratio_true,options)

    currentFolder = fileparts(mfilename('fullpath'));
    testSubFolder = fullfile(currentFolder, '..', 'methods');
    addpath(testSubFolder);

    folder_path = '../data_TrueSol';
    file_pattern = sprintf('type=%s_ratio=%.2f_Nx=%d_Nt=*_d=%.3f_lam1=%.2f_lam2=%.2f.mat',...
        problem_paras.type, mesh_ratio_true, Nx_true, problem_paras.d1,...
        problem_paras.lam1, problem_paras.lam2);
    
    files = dir(fullfile(folder_path, file_pattern));
    
    if isempty(files)
        fprintf('No matching files found.\n');
    else
        truedata = fullfile(folder_path, files(1).name);
        fprintf('Loading file: %s\n', truedata);
        load(truedata);
    end

    mult = Nx_true/Nx;

    True_piece1 = zeros(1,Nx+1);
    True_piece2 = zeros(1,Nx+1);
    
    for i = 0:Nx
        True_piece1(i+1) = TrueSol_last_row1(1+i*mult);
        True_piece2(i+1) = TrueSol_last_row2(1+i*mult);
    end

    if nargin < 6
        % Time of PCM
        tic;
        [last_row_pcm1,last_row_pcm2,~,~,~,~,L] = FDPCM(Nx,Nt,problem_paras);
        time_pcm = toc;
        % Err of PCM
        intF1 = last_row_pcm1 - True_piece1;
        intF2 = last_row_pcm2 - True_piece2;
        err_pcm = sqrt((intF1*intF1'+intF2*intF2')*2*L/Nx);
        % Time of RBF
        tic;
        [last_row_rbf1,last_row_rbf2,~,~,~,~,~] = RBCM(Nx,Nt,problem_paras);
        time_rbf = toc;
        % Err of PCM
        intF1 = last_row_rbf1 - True_piece1;
        intF2 = last_row_rbf2 - True_piece2;
        err_rbf = sqrt((intF1*intF1'+intF2*intF2')*2*L/Nx);
    else
        % Time of PCM
        tic;
        [last_row_pcm1,last_row_pcm2,~,~,~,~,L] = FDPCM(Nx,Nt,problem_paras,options);
        time_pcm = toc;
        % Err of PCM
        intF1 = last_row_pcm1 - True_piece1;
        intF2 = last_row_pcm2 - True_piece2;
        err_pcm = sqrt((intF1*intF1'+intF2*intF2')*2*L/Nx);
        % Time of RBF
        tic;
        [last_row_rbf1,last_row_rbf2,~,~,~,~,~] = RBCM(Nx,Nt,problem_paras,options);
        time_rbf = toc;
        % Err of PCM
        intF1 = last_row_rbf1 - True_piece1;
        intF2 = last_row_rbf2 - True_piece2;
        err_rbf = sqrt((intF1*intF1'+intF2*intF2')*2*L/Nx);
    end

    % Print the results
    func1_name = 'FDPCM';
    func2_name = 'RBCM';
    
    % Print the results with aligned columns
    format_str = 'The error and time of %-*s are: %10.6f, %10.6f\n';
    max_func_name_length = max(length(func1_name), length(func2_name));
    
    fprintf(format_str, max_func_name_length, func1_name, err_pcm, time_pcm);
    fprintf(format_str, max_func_name_length, func2_name, err_rbf, time_rbf);
end
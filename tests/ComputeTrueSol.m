function ComputeTrueSol(Nx,mesh_ratio,problem_paras,options)
    
    currentFolder = fileparts(mfilename('fullpath'));
    testSubFolder = fullfile(currentFolder, '..', 'methods');
    addpath(testSubFolder);
    TruncationTechFolder = fullfile(currentFolder, '..', 'methods/FDPCM_subfunc');
    addpath(TruncationTechFolder);

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

    Nt = round(problem_paras.T / ((2*L/Nx)^2/mesh_ratio));

    if nargin < 4
        [TrueSol_last_row1,TrueSol_last_row2,~,~,~,~,~] = FDPCM(Nx,Nt,problem_paras);
    else
        [TrueSol_last_row1,TrueSol_last_row2,~,~,~,~,~] = FDPCM(Nx,Nt,problem_paras,options);
    end
    
    % Construct the desired destination folder path
    destination_folder = fullfile(currentFolder, '..', 'data_TrueSol');
    filename = sprintf('type=%s_ratio=%.2f_Nx=%d_Nt=%d_d=%.3f_lam1=%.2f_lam2=%.2f.mat',...
        problem_paras.type, mesh_ratio, Nx, Nt, problem_paras.d1, problem_paras.lam1, problem_paras.lam2);
    full_path = fullfile(destination_folder, filename);
    save(full_path,'TrueSol_last_row1','TrueSol_last_row2')
end
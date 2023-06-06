function test
    currentFolder = fileparts(mfilename('fullpath'));
    testSubFolder = fullfile(currentFolder, '..', 'methods');
    addpath(testSubFolder);

    problem_paras.T = 1;
    problem_paras.K = 1;
    problem_paras.sig1 = 0.8;
    problem_paras.sig2 = 0.3;
    problem_paras.r1 = 0.05;
    problem_paras.r2 = 0.05;
    problem_paras.d1 = 0.025;
    problem_paras.d2 = 0.025;
    problem_paras.lam1 = 0.25;
    problem_paras.lam2 = 0.2;
    problem_paras.A = [-2 2; 3 -3];
    

    problem_paras.type = "Merton";


%     Plot2D_compare(256,500,problem_paras)

    Plot3D("FDPCM",64,1000,problem_paras)



%     options.ifprint = "yes";

%     TruncationTechFolder = fullfile(currentFolder, '..', 'methods/FDPCM_subfunc');
%     addpath(TruncationTechFolder);
%     switch problem_paras.type
%         case {"Merton"}
%             L = TruncationTech(problem_paras.K,0,problem_paras.sig1,...
%                 problem_paras.sig2,problem_paras.r1,problem_paras.r2,...
%                 problem_paras.d1,problem_paras.d2,problem_paras.lam1,...
%                 problem_paras.lam2,problem_paras.A);
%         case {"Kou"}
%             L = TruncationTech(problem_paras.K,1,problem_paras.sig1,...
%                 problem_paras.sig2,problem_paras.r1,problem_paras.r2,...
%                 problem_paras.d1,problem_paras.d2,problem_paras.lam1,...
%                 problem_paras.lam2,problem_paras.A);
%     end
%     Nx = 256;
%     Nt = round(problem_paras.T / ((2*L/Nx)^2/0.8))
%     Compare_ErrTime(Nx,Nt,problem_paras,1024,0.8)



%     TruncationTechFolder = fullfile(currentFolder, '..', 'methods/FDPCM_subfunc');
%     addpath(TruncationTechFolder);
%     switch problem_paras.type
%         case {"Merton"}
%             L = TruncationTech(problem_paras.K,0,problem_paras.sig1,...
%                 problem_paras.sig2,problem_paras.r1,problem_paras.r2,...
%                 problem_paras.d1,problem_paras.d2,problem_paras.lam1,...
%                 problem_paras.lam2,problem_paras.A);
%         case {"Kou"}
%             L = TruncationTech(problem_paras.K,1,problem_paras.sig1,...
%                 problem_paras.sig2,problem_paras.r1,problem_paras.r2,...
%                 problem_paras.d1,problem_paras.d2,problem_paras.lam1,...
%                 problem_paras.lam2,problem_paras.A);
%     end
%     Kou_Err_Time_vec = zeros(4,15);
%     Merton_Err_Time_vec = zeros(4,15);
%     Mesh_vec = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6];
%     Nx = 256;
% 
%     for i = 1:15
%         mesh_ratio_true = Mesh_vec(i);
%         Nt = round(problem_paras.T / ((2*L/Nx)^2/mesh_ratio_true));
%         [Kou_Err_Time_vec(1,i), Kou_Err_Time_vec(3,i), Kou_Err_Time_vec(2,i),...
%             Kou_Err_Time_vec(4,i)] = Compare_ErrTime(Nx,Nt,problem_paras,1024,mesh_ratio_true);
%     end


    


%     Nx = 512;
%     mesh_ratio = 5;
%     [full_path, filename] = ComputeTrueSol(Nx, mesh_ratio, problem_paras);
%     fprintf('full_path: %s\n', full_path);
%     fprintf('filename: %s\n', filename);


%     fixwhich = "spatial";
%     ErrorOrder(fixwhich,4,problem_paras,1024,0.8)
    

%     Plot_Freeboundaries(256,1000,problem_paras)

end
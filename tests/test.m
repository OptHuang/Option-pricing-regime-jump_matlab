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


%     options.ifprint = "yes";
%     Nt = 1000;
%     Nx = 32;
%     Compare_ErrTime(Nx,Nt,problem_paras,64,5)


    Nx = 512;
    mesh_ratio = 5;
    [full_path, filename] = ComputeTrueSol(Nx, mesh_ratio, problem_paras);
    fprintf('full_path: %s\n', full_path);
    fprintf('filename: %s\n', filename);


%     fixwhich = "ratio";
%     options.pcm_eps = 1e-7;
%     ErrorOrder(fixwhich,3,problem_paras,128,5,options)

end
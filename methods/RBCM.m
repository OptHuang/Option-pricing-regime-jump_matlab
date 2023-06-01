function [last_row1,last_row2,full_matrix1,full_matrix2,fb1,fb2,L] = RBCM(Nx,...
    Nt,problem_paras,options)
    % This function uses the method FDPCM to solve the American option
    % pricing under regime-switching with jumps.

    % last_row1: the value of option price under regime 1 at expiry time T,
    % and it is also the last row of full_matrix1.
    %
    % last_row2: the value of option price under regime 2 at expiry time T,
    % and it is also the last row of full_matrix2.
    %
    % full_matrix1: the solution matrix (size: (Nt+1)*(Nx+1)) of option
    % price under regime 1.
    %
    % full_matrix2: the solution matrix (size: (Nt+1)*(Nx+1)) of option
    % price under regime 2.
    %
    % Nx: the block number of the spatial partition.
    %
    % Nt: the block number of the temporal partition.
    %
    % problem_paras: structure of parameters for problem setup, including
    % "problem_paras".
    %
    % options: structure of parameters for PCM algorithm, including
    % "options.pcm_eps", "options.pcm_nu", "options.pcm_mu",
    % "options.pcm_vrho"; if we do not pass the "options", it will use the 
    % default settings.
    
    % Add the path "FDPCM_subfunc"
    currentFolder = fileparts(mfilename("fullpath"));
    rbcmSubFolder = fullfile(currentFolder, "RBCM_subfunc");
    addpath(rbcmSubFolder);

    % Set the option of the "options" to default values if it is not
    % supplied.
    if nargin < 4
        options = struct();
    end
    if isfield(options,"rbf_eps")
        rbf_eps = options.rbf_eps;
    else
        rbf_eps = 1e-6;
    end

    % Set parameters of the problem
    T = problem_paras.T;
    K = problem_paras.K;
    sig1 = problem_paras.sig1;
    sig2 = problem_paras.sig2;
    r1 = problem_paras.r1;
    r2 = problem_paras.r2;
    d1 = problem_paras.d1;
    d2 = problem_paras.d2;
    lam1 = problem_paras.lam1;
    lam2 = problem_paras.lam2;
    A = problem_paras.A;
    type = problem_paras.type;
    dt = T/Nt;

    switch type
        case {"Merton"}
            L = TruncationTech(K,0,sig1,sig2,r1,r2,d1,d2,lam1,lam2,A);
        case {"Kou"}
            L = TruncationTech(K,1,sig1,sig2,r1,r2,d1,d2,lam1,lam2,A);
    end

    [mu,delt,kappaM,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2] = Mertonpara(sig1,...
                sig2,r1,r2,d1,d2,lam1,lam2,A);
    [p,q,a1,a2,kappaK,nK1,nK2,mK1,mK2] = Koupara(sig1,sig2,...
                r1,r2,d1,d2,lam1,lam2,A);
    switch type
        case {"Merton"}
            [Coefmatrix,A1] = RBFassemblemat(1,0,dt,Nx,A,sig1,sig2,r1,r2,...
                d1,d2,lam1,lam2,kappaM,kappaK,mu,delt,p,q,a1,a2,L,rbf_eps);
        case {"Kou"}
            [Coefmatrix,A1] = RBFassemblemat(1,1,dt,Nx,A,sig1,sig2,r1,r2,...
                d1,d2,lam1,lam2,kappaM,kappaK,mu,delt,p,q,a1,a2,L,rbf_eps);
    end

    [Lmat,Umat] = lu(Coefmatrix);
    [full_matrix1,full_matrix2] = SolveRBF(Lmat,Umat,A1,L,Nx,Nt,K);
    
    last_row1 = full_matrix1(end,:);
    last_row2 = full_matrix2(end,:);

    [~,~,fb1,fb2] = FreeBoundry(full_matrix1,full_matrix2,K,L,Nx,Nt);

    rmpath(rbcmSubFolder)
end
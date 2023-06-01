function [last_row1,last_row2,full_matrix1,full_matrix2,fb1,fb2,L] = FDPCM(Nx,...
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
    % "options.pcm_vrho", "options.ifprint"; if we do not pass the
    % "options", it will use the default settings.


    % Add the path "FDPCM_subfunc"
    currentFolder = fileparts(mfilename("fullpath"));
    fdpcmSubFolder = fullfile(currentFolder, "FDPCM_subfunc");
    addpath(fdpcmSubFolder);

    % Set the option of the "options" to default values if it is not
    % supplied.
    if nargin < 4
        options = struct();
    end
    if isfield(options,"pcm_eps")
        pcm_eps = options.pcm_eps;
    else
        pcm_eps = get_default("pcm_eps");
    end
    if isfield(options,"pcm_nu")
        pcm_nu = options.pcm_nu;
    else
        pcm_nu = get_default("pcm_nu");
    end
    if isfield(options,"pcm_mu")
        pcm_mu = options.pcm_mu;
    else
        pcm_mu = get_default("pcm_mu");
    end
    if isfield(options,"pcm_vrho")
        pcm_vrho = options.pcm_vrho;
    else
        pcm_vrho = get_default("pcm_vrho");
    end
    if isfield(options,"ifprint") 
        % To decide whether to print the iterator in each iteration (to
        % check the procedure); if ifprint="yes", then print; if
        % ifprint="no", then not print.
        ifprint = options.ifprint;
    else
        ifprint = get_default("ifprint");
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

    switch type
        case {"Merton"}
            L = TruncationTech(K,0,sig1,sig2,r1,r2,d1,d2,lam1,lam2,A);
        case {"Kou"}
            L = TruncationTech(K,1,sig1,sig2,r1,r2,d1,d2,lam1,lam2,A);
    end

    % Initialize.
    full_matrix1 = zeros(Nt+1,Nx+1);
    full_matrix2 = zeros(Nt+1,Nx+1);
    x = linspace(-L,L,Nx+1);
    S = exp(x);
    dx = 2*L/Nx;
    dt = T/Nt;
    x0 = -L;
    t0 = 0;

    % Assemble matrices and solve the problem
    if ifprint=="yes"
        switch type
            case {"Merton"}
                [mu,delt,~,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2] = Mertonpara(sig1,...
                    sig2,r1,r2,d1,d2,lam1,lam2,A);
                w1_star = exp(-nM1*x).*max(K-exp(x),0);
                Phi1_star = w1_star(2:1:Nx).';
                w2_star = exp(-nM2*x).*max(K-exp(x),0);
                Phi2_star = w2_star(2:1:Nx).';
                Phi_pre = [Phi1_star;Phi2_star];
                [Merton1,Merton2] = AssembleRMerton(Nx,delt,mu_star1,mu_star2,dx,C1,C2);
                [B, d1t0, d2t0, F0t0, Phi_star_t0] = AssembleBasicMatMerton(dt,...
                    dx,Nx,sig1,sig2,nM1,nM2,mM1,mM2,A,lam1,lam2,L,K,...
                    mu_star1,mu_star2,mu,delt);
    
                for i = 2:Nt+1
                    fprintf("%d\n", i)
                    [F, Phi_star] = AssemblecurrMerton(i-1,dt,Nx,mM1,mM2,lam1,...
                        lam2,Phi_pre,Merton1,Merton2,B,d1t0,d2t0,F0t0,Phi_star_t0);
                    len = size(F,1)/2;
                    Phi_cur = zeros(len*2,1);
                    Phi_start = Phi_pre-Phi_star;
                    for j = 1:2
                        F_temp = F((j-1)*len+1:j*len);
                        B_temp = B((j-1)*len+1:j*len,(j-1)*len+1:j*len);
                        Phi_start_temp = Phi_start((j-1)*len+1:j*len);
                        Phi_cur((j-1)*len+1:j*len) = PCMPro(pcm_eps,pcm_nu,pcm_mu,...
                            pcm_vrho,B_temp,F_temp,Phi_start_temp);
                    end
                    full_matrix1(i,2:Nx) = Phi_cur(1:Nx-1).'+Phi_star(1:Nx-1).';
                    full_matrix2(i,2:Nx) = Phi_cur(Nx:2*Nx-2).'+Phi_star(Nx:2*Nx-2).';
                    [V1,V2] = RevValueMerton(i-2,Nx,Phi_cur,Phi_star,mM1,mM2,nM1,nM2,t0,x0,dt,dx);
                    full_matrix1(i,2:Nx) = V1.';
                    full_matrix2(i,2:Nx) = V2.';
                    Phi_pre = Phi_cur+Phi_star;
                end
            case {"Kou"}
                [p,q,a1,a2,~,nK1,nK2,mK1,mK2] = Koupara(sig1,sig2,...
                    r1,r2,d1,d2,lam1,lam2,A);
                w1_star = exp(-nK1*x).*max(K-exp(x),0);
                Phi1_star = w1_star(2:1:Nx).';
                w2_star = exp(-nK2*x).*max(K-exp(x),0);
                Phi2_star = w2_star(2:1:Nx).';
                Phi_pre = [Phi1_star;Phi2_star];
                [Kou1,Kou2] = AssembleRKou(Nx,nK1,nK2,dx,p,q,a1,a2);
                [B, d1t0, d2t0, F0t0, Phi_star_t0] = AssembleBasicMatKou(dt,...
                    dx,Nx,sig1,sig2,nK1,nK2,mK1,mK2,A,lam1,lam2,L,K,p,q,a1,a2);
    
                for i = 2:Nt+1
                    fprintf("%d\n", i)
                    [F, Phi_star] = AssemblecurrKou(i-1,dt,Nx,mK1,mK2,lam1,...
                        lam2,Phi_pre,Kou1,Kou2,B, d1t0, d2t0, F0t0, Phi_star_t0);
                    len = size(F,1)/2;
                    Phi_cur = zeros(len*2,1);
                    Phi_start = Phi_pre-Phi_star;
                    for j = 1:2
                        F_temp = F((j-1)*len+1:j*len);
                        B_temp = B((j-1)*len+1:j*len,(j-1)*len+1:j*len);
                        Phi_start_temp = Phi_start((j-1)*len+1:j*len);
                        Phi_cur((j-1)*len+1:j*len) = PCMPro(pcm_eps,pcm_nu,pcm_mu,...
                            pcm_vrho,B_temp,F_temp,Phi_start_temp);
                    end
                    full_matrix1(i,2:Nx) = Phi_cur(1:Nx-1).'+Phi_star(1:Nx-1).';
                    full_matrix2(i,2:Nx) = Phi_cur(Nx:2*Nx-2).'+Phi_star(Nx:2*Nx-2).';
                    [V1,V2] = RevValueKou(i-2,Nx,Phi_cur,Phi_star,mK1,mK2,nK1,nK2,t0,x0,dt,dx);
                    full_matrix1(i,2:Nx) = V1.';
                    full_matrix2(i,2:Nx) = V2.';
                    Phi_pre = Phi_cur+Phi_star;
                end
        end
    else
        switch type
            case {"Merton"}
                [mu,delt,~,nM1,nM2,mM1,mM2,mu_star1,mu_star2,C1,C2] = Mertonpara(sig1,...
                    sig2,r1,r2,d1,d2,lam1,lam2,A);
                w1_star = exp(-nM1*x).*max(K-exp(x),0);
                Phi1_star = w1_star(2:1:Nx).';
                w2_star = exp(-nM2*x).*max(K-exp(x),0);
                Phi2_star = w2_star(2:1:Nx).';
                Phi_pre = [Phi1_star;Phi2_star];
                [Merton1,Merton2] = AssembleRMerton(Nx,delt,mu_star1,mu_star2,dx,C1,C2);
                [B, d1t0, d2t0, F0t0, Phi_star_t0] = AssembleBasicMatMerton(dt,...
                    dx,Nx,sig1,sig2,nM1,nM2,mM1,mM2,A,lam1,lam2,L,K,...
                    mu_star1,mu_star2,mu,delt);
    
                for i = 2:Nt+1
                    [F, Phi_star] = AssemblecurrMerton(i-1,dt,Nx,mM1,mM2,lam1,...
                        lam2,Phi_pre,Merton1,Merton2,B,d1t0,d2t0,F0t0,Phi_star_t0);
                    len = size(F,1)/2;
                    Phi_cur = zeros(len*2,1);
                    Phi_start = Phi_pre-Phi_star;
                    for j = 1:2
                        F_temp = F((j-1)*len+1:j*len);
                        B_temp = B((j-1)*len+1:j*len,(j-1)*len+1:j*len);
                        Phi_start_temp = Phi_start((j-1)*len+1:j*len);
                        Phi_cur((j-1)*len+1:j*len) = PCMPro(pcm_eps,pcm_nu,pcm_mu,...
                            pcm_vrho,B_temp,F_temp,Phi_start_temp);
                    end
                    full_matrix1(i,2:Nx) = Phi_cur(1:Nx-1).'+Phi_star(1:Nx-1).';
                    full_matrix2(i,2:Nx) = Phi_cur(Nx:2*Nx-2).'+Phi_star(Nx:2*Nx-2).';
                    [V1,V2] = RevValueMerton(i-2,Nx,Phi_cur,Phi_star,mM1,mM2,nM1,nM2,t0,x0,dt,dx);
                    full_matrix1(i,2:Nx) = V1.';
                    full_matrix2(i,2:Nx) = V2.';
                    Phi_pre = Phi_cur+Phi_star;
                end
            case {"Kou"}
                [p,q,a1,a2,~,nK1,nK2,mK1,mK2] = Koupara(sig1,sig2,...
                    r1,r2,d1,d2,lam1,lam2,A);
                w1_star = exp(-nK1*x).*max(K-exp(x),0);
                Phi1_star = w1_star(2:1:Nx).';
                w2_star = exp(-nK2*x).*max(K-exp(x),0);
                Phi2_star = w2_star(2:1:Nx).';
                Phi_pre = [Phi1_star;Phi2_star];
                [Kou1,Kou2] = AssembleRKou(Nx,nK1,nK2,dx,p,q,a1,a2);
                [B, d1t0, d2t0, F0t0, Phi_star_t0] = AssembleBasicMatKou(dt,...
                    dx,Nx,sig1,sig2,nK1,nK2,mK1,mK2,A,lam1,lam2,L,K,p,q,a1,a2);
    
                for i = 2:Nt+1
                    [F, Phi_star] = AssemblecurrKou(i-1,dt,Nx,mK1,mK2,lam1,...
                        lam2,Phi_pre,Kou1,Kou2,B, d1t0, d2t0, F0t0, Phi_star_t0);
                    len = size(F,1)/2;
                    Phi_cur = zeros(len*2,1);
                    Phi_start = Phi_pre-Phi_star;
                    for j = 1:2
                        F_temp = F((j-1)*len+1:j*len);
                        B_temp = B((j-1)*len+1:j*len,(j-1)*len+1:j*len);
                        Phi_start_temp = Phi_start((j-1)*len+1:j*len);
                        Phi_cur((j-1)*len+1:j*len) = PCMPro(pcm_eps,pcm_nu,pcm_mu,...
                            pcm_vrho,B_temp,F_temp,Phi_start_temp);
                    end
                    full_matrix1(i,2:Nx) = Phi_cur(1:Nx-1).'+Phi_star(1:Nx-1).';
                    full_matrix2(i,2:Nx) = Phi_cur(Nx:2*Nx-2).'+Phi_star(Nx:2*Nx-2).';
                    [V1,V2] = RevValueKou(i-2,Nx,Phi_cur,Phi_star,mK1,mK2,nK1,nK2,t0,x0,dt,dx);
                    full_matrix1(i,2:Nx) = V1.';
                    full_matrix2(i,2:Nx) = V2.';
                    Phi_pre = Phi_cur+Phi_star;
                end
        end
    end
    

    for i = 2:Nx
        full_matrix1(1,i) = V_star(K,S(i));
        full_matrix2(1,i) = V_star(K,S(i));
    end
    % 右侧
    for i = 1:Nt+1
        full_matrix1(i,Nx+1) = 0;
        full_matrix2(i,Nx+1) = 0;
    end
    % 左侧
    for i = 1:Nt+1
        full_matrix1(i,1) = V_star(K,exp(-L));
        full_matrix2(i,1) = V_star(K,exp(-L));
    end

    last_row1 = full_matrix1(end,:);
    last_row2 = full_matrix2(end,:);

    [~,~,fb1,fb2] = FreeBoundry(full_matrix1,full_matrix2,K,L,Nx,Nt);

    rmpath(fdpcmSubFolder);
end
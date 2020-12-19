function result =  run_pgplvm_on_gaus(xxtrue,yytrue, ctrs)
    
    nf = size(xxtrue,2); % number of latent dimensions
    xplds = run_plds(yytrue,nf)';
    xpldsmat = align_xtrue(xplds,xxtrue); % align the estimate with the true latent variable.
    xinit = xplds;
    xinitmat = xpldsmat;
    ff = ffun_bump(xxtrue, ctrs);
    
    setopt.sepx_flag = 0;
    setopt.sigma2_init = 3;
    setopt.sigma2_end = min([0.1,setopt.sigma2_init]); % initial noise variance
    setopt.lr = 0.95; % learning rate
    setopt.latentTYPE = 2; % kernel for the latent, 1. AR1, 2. SE
    setopt.ffTYPE = 2; % kernel for the tuning curve, 1. AR1, 2. SE, 3. linear,4. SE_len
    setopt.initTYPE = 2; % initialize latent: 1. use PLDS init; 2. use random init; 3. true xxtrue
    setopt.la_flag = 1; % 1. no la; 2. standard la; 3. decoupled la,  obsolete
    setopt.rhoxx = 1; % rho for Kxx
    setopt.lenxx = 5; % len for Kxx
    setopt.rhoff = 1; % rho for Kff
    setopt.lenff = 5; % median(vec(range(xinit))); % len for Kff
    setopt.lenff_ratio = 1; % len ratio for Kff
    setopt.b = 0; % obsolete
    setopt.r = 1; % obsolete
    setopt.nsevar = 1; % obsolete
    setopt.hypid = [2,3,4]; % 1. rho for Kxx; 2. len for Kxx; 3. rho for Kff; 4. len for Kff; 5. sigma2 (annealing it instead of optimizing it)
    setopt.xinit = xinit; % for initialization purpose
    setopt.niter = 500; % number of iterations
    setopt.opthyp_flag = 0;
    setopt.xplds = xplds;
    setopt.xpldsmat = xpldsmat;
    
    result = pgplvm_la(yytrue, xxtrue, ff, setopt);

    
    
end


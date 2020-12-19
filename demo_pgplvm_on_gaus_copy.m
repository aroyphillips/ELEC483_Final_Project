function result = demo_pgplvm_on_gaus(xx,yytrue,ff)
%
% Adapted from Anqi Wu's tutorial script demo3_2DBump_phil illustrating P-GPLVM for 2-dimensional latent variable
% with tuning curves generated from 2D Gaussian bumps.
% Changes include setting experiment pameters (nneur, etc) making latent kernel SE, and removing outdated
% function calls with missing dependencies (eg covSEiso_sp, which leads to missing function "realbasis_nD_curv_kron".

set(0,'DefaultFigureVisible','off'); % added to stop taking focus
% Initialize paths
initpaths;

% dataSource = 1; %0 for matlab generated, 1 for python generated
% 
% % Load data
% 
% if dataSource == 0
%     datasetname = 'simdatadir/simdata3_phillips.mat';  % name of dataset
%     if ~exist(datasetname,'file') % Create simulated dataset if necessary
%         fprintf('Creating simulated dataset: ''%s''\n', datasetname);
%         mkSimData3_2DBump;
%     end
%     load(datasetname);
%     xxtrue = simdata.latentVariable;
%     yytrue = simdata.spikes;
%     fftrue = simdata.spikeRates;
% elseif dataSource ==1
%     xxtrue = readmatrix('latent_bump.csv')'; %simdata.latentVariable;
%     yytrue = readmatrix('sim_spikes_bump.csv')'; %simdata.spikes;
%     fftrue = readmatrix('logFR_bump.csv'); %simdata.spikeRates;
% else
%     warning("dataSource must be set to 0 for matlab generated or 1 for python generated");
% end
%     

% Get sizes and spike counts
[nt0,nneur] = size(yytrue); % nt: number of time points; nneur: number of neurons
nf = size(xx,2); % number of latent dimensions

%% == 1. Compute baseline estimates ====

% Initialize the log of spike rates with the square root of spike counts.
ffmat = sqrt(yytrue);

% % Compute LLE
% xlle = lle(ffmat',20,nf)';
% xllemat = align_xtrue(xlle,xx); % align the estimate with the true latent variable.

% Compute PPCA
% options = fgplvmOptxions('ftc');
xppca = pca(ffmat,nf); % this is dead for some reason
xppcamat = align_xtrue(xppca,xx); % align the estimate with the true latent variable.

% Compute Poisson Linear Dynamic System (PLDS)
xplds = run_plds(yytrue,nf)';
xpldsmat = align_xtrue(xplds,xx); % align the estimate with the true latent variable.

% truncate into segments, re-organize the data
if nt0>=2000
    nt = find_nt(nt0);
else
    nt = nt0;
end
ntr = nt0/nt;
xinit = xplds;
yy_all = permute(reshape(yytrue,nt,ntr,1,[]),[2,3,1,4]);
xinit_all = permute(reshape(xinit,nt,ntr,1,[]),[2,3,1,4]);
sid = 1:ntr;
tid = 1;
yy = reshape(yy_all(sid,tid,:,:),[],nt,nneur);
xinit = reshape(xinit_all(sid,tid,:,:),[],nt,nf);
if sum(vec(abs(xinit)))==0
    xinit = randn(size(xinit))*1e-5;
end

seg_id = zeros(length(sid),length(tid));
for ii=1:length(sid)
    for jj=1:length(tid)
        seg_id(ii,jj,:) = tid(jj);
    end
end
seg_id = vec(seg_id);
%% == 2. Compute P-GPLVM ====

% Set up options
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
setopt.xinit = xplds; % for initialization purpose
setopt.niter = 50; % number of iterations
setopt.xplds = xplds;
setopt.xpldsmat = xpldsmat;
setopt.opthyp_flag = 0;


% Compute P-GPLVM with Laplace Approximation
result = pgplvm_la_updated(yy,xx,ff,setopt);


%% == 3. Plot tuning curves ====
figure(2)
xgrid = gen_grid([min(result_la.xxsampmat(:,1)) max(result_la.xxsampmat(:,1)); min(result_la.xxsampmat(:,2)) max(result_la.xxsampmat(:,2))],50,nf); % x grid for plotting tuning curves
fftc = exp(get_tc(result_la.xxsampmat,result_la.ffmat,xgrid,result_la.rhoff,result_la.lenff));
fftc_true = simdata.tuningCurve;

chosenNeurs = [1,4,14];
for ii=choseNeurs
    fg = reshape(fftc_true(:,ii),50,[]);
    xg = gen_grid([-6 6],50,1);
    subplot(131),cla,hold on
    contourf(xg,xg,fg)
    plot(xx(:,1),xx(:,2),'g.-')
    hold off
    title(['tuning curve for neuron ' num2str(ii)]);
    xlabel('x');
    set(gca,'fontsize',15)
    
    subplot(132),surf(xg,xg,reshape(fftc(:,ii),50,[])),title('estimated tc')
    set(gca,'fontsize',15)
    subplot(133),surf(xg,xg,reshape(fftc_true(:,ii),50,[])),title('true tc')
    set(gca,'fontsize',15)
    drawnow,pause
end

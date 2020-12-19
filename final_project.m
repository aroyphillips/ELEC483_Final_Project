%% Roy Phillips ELEC 483 Final Project:
%%% testing algorithm on simulated data
clear;clc;close all;

%% load simulated grid neuron data

niter = 500;
xtrueSin = readmatrix('latent_sin.csv')';
spikesSin = readmatrix('sim_spikes_sin.csv')';
sinParams = readmatrix('sin_params.csv');
omega = sinParams(:,2);
phi = sinParams(:,3);

kt= readmatrix('kt.csv');

%%
ff_sin = ffun_sin(xtrueSin, omega, phi);

result_sin = demo_pgplvm(xtrueSin, spikesSin, ff_sin, niter);

%% plot sin tuning tuning curve

chosenSinNeurs = [7,9,19];
colors = {'b','r','g'};
figure;
hold on
idx=1;
for neur = chosenSinNeurs
    hold on
    plot(result_sin.xxsampmat, exp(result_sin.ffmat(:,neur)), colors{idx});
    legend({'Neuron 6', 'Neuron 8', 'Neuron 18'})
    idx = idx+1;
end
title("estimated tuning curves for grid cells")
ylabel("firing rate")
xlabel("location")

hold off


%% load simulated gaussian bump data

% Load data
isDataSourcePython = 1;
[xxtrue, yytrue, fftrue, centers] = load_gaus_data(isDataSourcePython); % loads my simulated data from either matlab or python

result = demo_pgplvm(xxtrue, yytrue, fftrue, 5000); 




%% Plot bump tuning curves ====
figure;
xgrid = gen_grid([min(result.xxsampmat(:,1)) max(result.xxsampmat(:,1)); min(result.xxsampmat(:,2)) max(result.xxsampmat(:,2))],50,2); % x grid for plotting tuning curves
fftc = exp(get_tc(result.xxsampmat,result.ffmat,xgrid,result.rhoff,result.lenff));
ff = exp(get_tc(xxtrue,fftrue,xgrid,result.rhoff,result.lenff));

chosenNeurs = [1,14,20];
idx = 1;
for ii=chosenNeurs
      fg = reshape(ff(:,ii),50,[]);
      xg = gen_grid([-6 6],50,1);
      subplot(3,1,idx),pcolor(xg,xg,reshape(fftc(:,ii),50,[])),title(['estimated tuning curve for neuron' num2str(ii)])
      set(gca,'fontsize',15)
    drawnow
    idx = idx+1;
end

pause;
%% functions
function ff = ffun_sin(x, omega, phi)
    %%% sinusoidal latent function
    %%% x is size (T,P) matrix, omega is (N,P), phi is (N,P)
    %%% returns mapping ff = 3*sin(omega*x+phi) = ln firing rates, size=(T,N)
    ff = (3*sin(omega*x'+phi)-2)';
end
    
function [xxtrue, yytrue, fftrue, centers] = load_gaus_data(dataSource)
    % given a boolean (0=MATLAB) (1=Python), loads the corresponding
    % simulated data (both datasets simulated by Roy Phillips with the same centers, different latent).
    if dataSource == 0
        datasetname = 'simdatadir/simdata3_phillips.mat';  % name of dataset
        if ~exist(datasetname,'file') % Create simulated dataset if necessary
            fprintf('Creating simulated dataset: ''%s''\n', datasetname);
            mkSimData3_2DBump;
        end
        load(datasetname);
        xxtrue = simdata.latentVariable;
        yytrue = simdata.spikes;
        fftrue = simdata.spikeRates;
        centers = simdata.centers;
    elseif dataSource ==1
        xxtrue = readmatrix('latent_bump.csv')'; 
        yytrue = readmatrix('sim_spikes_bump.csv')';
        fftrue = readmatrix('logFR_bump.csv');
        centers = readmatrix('gaus_centers.csv');
    else
        warning("dataSource must be set to 0 for matlab generated or 1 for python generated");
    end
end

function ff_bump = ffun_bump(xx, ctrs, ell, rho)
    %% xx is (T,P), ctrs (N, P), ell & rho are hyperparametes
    if nargin <3
        % default vals
        ell = 1; % standard deviation
        rho = 20; % marginal variance
    end

    K = sq_dist(xx'/ell,ctrs'/ell);  % cross covariances Kxz

    ff_bump = rho*exp(-K/2);
        
end

function [xxtrue, yytrue, fftrue, centers] = load_sin_data()
    % given a boolean (0=MATLAB) (1=Python), loads the corresponding
    % simulated data (both datasets simulated by Roy Phillips with the same centers, different latent).
    if dataSource == 0
        datasetname = 'simdatadir/simdata3_phillips.mat';  % name of dataset
        if ~exist(datasetname,'file') % Create simulated dataset if necessary
            fprintf('Creating simulated dataset: ''%s''\n', datasetname);
            mkSimData3_2DBump;
        end
        load(datasetname);
        xxtrue = simdata.latentVariable;
        yytrue = simdata.spikes;
        fftrue = simdata.spikeRates;
        centers = simdata.centers;
    elseif dataSource ==1
        xxtrue = readmatrix('latent_bump.csv')'; 
        yytrue = readmatrix('sim_spikes_bump.csv')';
        fftrue = readmatrix('logFR_bump.csv');
        centers = readmatrix('gaus_centers.csv');
    else
        warning("dataSource must be set to 0 for matlab generated or 1 for python generated");
    end
end

    
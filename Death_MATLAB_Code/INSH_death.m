% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% MATLAB Script to find optimal Bayesian designs via  %
% the Induced Natural Selection Heuristic (INSH)      %
%                                                     %
% Author: David Price                                 %
% Date: 4th May, 2018                                 %
%                                                     %
% Details of the method can be found in:              %
%                                                     %
% An Induced Natural Selection Heuristic for Finding  %
%   Optimal Bayesian Experimental Designs (2018),     %
%   Computational Statistics and Data Analysis,       %
%   ():x-y                                            %
%                                                     %
%                                                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %

rng(1)

ntimes = 1;

number_of_data=50000;


% % % % % % % % % % % % % % % % % % % % % % % % %
% % % Define prior distribution parameters% % % %
% % % % % % % % % % % % % % % % % % % % % % % % %
prior_b1=-0.005;
% Prior std
sd_b1=sqrt(0.01);

upper_b1=1.5;
lower_b1=0.6;
Tmax=10;

% % % % % % % % % % % % % % % % % % % % % %
% % % Define experimental parameters% % % %
% % % % % % % % % % % % % % % % % % % % % %
% initial number of infected
N=50;

% To pass to the sampling function for new designs at each wave
stepsize = 0.05;

T = [stepsize, 10];

tol=0.25*ntimes;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % %  Define prior distribution across grid points % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Grid (alpha, rho)
n=100;

prior_grid=(lower_b1):((upper_b1-lower_b1)/n):upper_b1;

prior_values = lognrnd(prior_b1*ones(number_of_data,1),sd_b1);
prior_probabilities = histc(prior_values, prior_grid)/number_of_data;

if size(prior_probabilities,1)~=1
    prior_probabilities=prior_probabilities';
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % Set CE-ABCDE convergence parameters % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

switch ntimes
    case 1
        W=8;
        m=[3*ones(1,W/2), 5*ones(1,W/2)];
        ntokeep = [10*ones(1,W/2), 6*ones(1,W/2)];
        ninit=20;
    case 2
        W=10;
        m=[3*ones(1,W/2), 5*ones(1,W/2)];
        ntokeep = [20*ones(1,W/2), 12*ones(1,W/2)];
        ninit=50;
    case 3
        W=16;
        m=[3*ones(1,W/2), 5*ones(1,W/2)];
        ntokeep = [20*ones(1,W/2), 12*ones(1,W/2)];
        ninit=120;
    case 4
        W=20;
        m=[3*ones(1,W/2), 5*ones(1,W/2)];
        ntokeep = [20*ones(1,W/2), 12*ones(1,W/2)];
        ninit=250;
    case 6
        W=30;
        m=[3*ones(1,W/2), 5*ones(1,W/2)];
        ntokeep = [25*ones(1,W/2), 15*ones(1,W/2)];
        ninit=400;
    case 8
        W=50;
        m=[3*ones(1,W/2), 5*ones(1,W/2)];
        ntokeep = [25*ones(1,W/2), 15*ones(1,W/2)];
        ninit=600;
        tol = 1.5;
end


sample_sd= 0.1*ones(1,ntimes);

%%
% Initial wave of designs is these times

current_wave_designs = round(sort(unifrnd(0,10,[ninit,ntimes]),2)/stepsize)*stepsize;


tic;

keep_all = [];

for w=1:W
    
    sprintf('Wave %u.', w)
    % sprintf('No. current wave designs: %u', size(current_wave_designs,1))
    % sprintf('No. unique time points: %u', length(unique(current_wave_designs)))
    
    data = zeros(number_of_data,length(unique(current_wave_designs)));
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % Simulate data at each of the sparse observation times % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    parfor d=1:size(prior_values,1)
        data(d,:)=death_simulate(prior_values(d), unique(current_wave_designs)', N);
    end
    
    kld=zeros(size(current_wave_designs,1),1);
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % Evaluate utility at each design using ABCdE approach  % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    parfor i=1:size(current_wave_designs,1)
        tmp_data = data(:,arrayfun(@(x)find(unique(current_wave_designs)==x,1),current_wave_designs(i,:)));
        
        [uni_data, num] = howmanyunique(tmp_data);
        
        tmp_kld=zeros(size(uni_data,1),1);
        std_tmp=std(tmp_data);
        
        for k=1:size(uni_data,1)
            
            keep = prior_values(ABCDE_chunk(tmp_data, uni_data(k,:), std_tmp)<tol);
            
            % IF not using ABCDE_chunk, the following line does the same operation directly
            % keep = prior_values(sum(bsxfun(@rdivide, abs(bsxfun(@minus, tmp_data, uni_data(k,:))), std_tmp),2) < tol);
            
            density = histc(keep, prior_grid);
            
            density=bsxfun(@rdivide, density, sum(density));
            
            if size(density,1)~=1
                density=density';
            end
            
            tmp_kld(k) = sum(nansum(log(density./prior_probabilities).*density))/length(density);
        end
        
        kld(i) = dot(tmp_kld, num/sum(num));
        
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % Store all sampled design points and corresponding utility % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    keep_all = [keep_all; [current_wave_designs, kld, w*ones(length(kld),1)]];
    
    
    sortutils= sort(kld, 'descend');
    util_cutoff = sortutils(ntokeep(w));
    update_wave = current_wave_designs(kld >= util_cutoff, :);
    
    % sprintf('Size update wave: %u',size(update_wave,1))
    
    % Include the optimal design thus far in the update wave so that we
    % keep exploring the space around the optimal, rather than hitting it
    % then moving away from it
    if(keep_all(keep_all(:,ntimes+1)==max(keep_all(:,ntimes+1)),end)~=w)
        update_wave = [update_wave; keep_all((keep_all(:,ntimes+1)==max(keep_all(:,ntimes+1))),1:ntimes)];
    end
    %  sprintf('Size update wave: %u',size(update_wave,1))
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % Sample new designs around accepted designs % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    current_wave_designs = [];
    for k=1:size(update_wave,1)
        current_wave_designs = [current_wave_designs; discrete_trunc_sample(update_wave(k,:), T, sample_sd, stepsize, m(w));];
    end
    
    % Add OED back in to be re-evaluated at next iteration
    current_wave_designs = [current_wave_designs; keep_all((keep_all(:,ntimes+1)==max(keep_all(:,ntimes+1))),1:ntimes)];
    
end

runtime = toc;
%



r=find(keep_all(:,ntimes+1)==max(keep_all(:,ntimes+1)));
opt_times = keep_all(r,1:ntimes)
opt_util = keep_all(r, ntimes+1);
v = unique(keep_all(:,end)); c = hist(keep_all(:,end),v);
ndesigns_per_wave = [v'; c]
total_designs_evaluated = size(keep_all,1);
%


sprintf('Run time: %f', runtime)

sprintf('No. designs considered: %u', total_designs_evaluated)

sprintf('Optimal occured in wave %u', unique(keep_all(find(keep_all(:,ntimes+1)==max(keep_all(:,ntimes+1))),end)))


filename = sprintf('INSH_death_t%d.mat', ntimes)
save(filename, 'keep_all', 'runtime','opt_times','total_designs_evaluated')

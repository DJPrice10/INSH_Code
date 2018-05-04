function new_design = discrete_trunc_sample(T, minmaxT, sd_T, stepsize, m)
% 
% sd_T = [0.1, 0.1, 0.1];
% T = [4.9 5 8];
% stepsize = 0.1;
% minmaxT = [0.1 10];
% m=10;
% 
% 
% % Pre-define cell array with one element for each obs time
% gridargs = cell(1,length(T));
% 
% % For each obs time, create a vector of times ± 10 sd's either direction
% % (i.e., create a boT/cube/etc. around each current obs time)
% for i=1:length(T)
%     gridargs{i} = (max(T(i)-8*sd_T(i), min(minmaxT)):stepsize:min(T(i)+8*sd_T(i), max(minmaxT)));
% end
% 
% % Pre-define a cell array of the outputs from ndgrid (the gridded times)
% outs = cell(1,length(T));
% 
% [outs{:}] = ndgrid(gridargs{:});
% 
% Ts = zeros(numel(outs{1}), length(T));
% 
% % Convert the array of times into a matrix with m rows (no. combinations)
% % and n columns (no. observations)
% for i=1:length(outs)
%     Ts(:,i) = outs{i}(:);
% end
% 
% % Get differences in each column
% tmp = diff(Ts, 1,2);
% 
% % Get rows that all times satisfy t1<t2<... constraint
% idx = ~any(tmp<=0, 2);
% 
% % Get those rows
% Ts = Ts(idx,:);
% 
% % Get probabilities from MV normal centred at current time, with sd_t, and
% % evaluate the probabilities on the grid specified in Ts
% probs = mvnpdf(Ts, T, sd_T);
% 
% % Sample a new design on the grid with weights given by probs
% new_design = Ts(randsample(1:size(Ts,1), m, 'true', probs), :);
% 
% % figure(2);
% % [~, density, X, Y] = kde2d(new_design(:,1:2), 2^6, [0.1 0.1], [10, 10]);
% % surf(X, Y, density)
% % view(0,90)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_design = -1*ones(m,length(T));
for i=1:m
    while( (any(any(new_design(i,:)<min(minmaxT))) || any(any(diff(new_design(i,:),1,2)<=0)) ) )
        new_design(i,:) = round( mvnrnd(T, sd_T, 1) /stepsize) * stepsize;
    end
end

function p = prt_friedman(PRT, models)
% Function to perform a Friedman test
%
% Inputs:
% -------
% PRT:     PRT struct
% models:  model IDs (int vector)
%
% Outputs:
% --------
% p:       p-value struct
%
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J.M. Monteiro


% Make initial checks
% -------------------------------------------------------------------------
if length(models) < 2
    error('Number of models must be greater than 1.');
end

for m = models
    if ~strcmp(PRT.model(models(1)).input.type, PRT.model(m).input.type)
        error('Models must be of the same type (classification or regression).');
    end
    if length(PRT.model(models(1)).output.fold) ~= length(PRT.model(m).output.fold)
        error('Models must have the same number of folds.');
    end
end


% Get stats for all folds in both models
% -------------------------------------------------------------------------
stats = cell(length(PRT.model(models(1)).output.fold), length(models));
i = 1;
for m = models
    for f = 1:length(PRT.model(m).output.fold)
        stats{f,i} =  PRT.model(m).output.fold(f).stats;
    end
    i = i+1;
end

% Friedman Test
% -------------------------------------------------------------------------
if strcmp(PRT.model(models(1)).input.type, 'classification')
    
    acc = nan(size(stats));
    b_acc = nan(size(stats));
    for i = 1:size(stats,1)
        for j = 1:size(stats,2)
            acc(i,j) = stats{i,j}.acc;
            b_acc(i,j) = stats{i,j}.b_acc;
        end
    end
    
    p.acc = friedman(acc, 1, 'off');
    p.b_acc = friedman(b_acc, 1, 'off');
    
elseif strcmp(PRT.model(models(1)).input.type, 'regression')
    
    corr = nan(size(stats));
    r2 = nan(size(stats));
    mse = nan(size(stats));
    nmse = nan(size(stats));
        
    for i = 1:size(stats,1)
        for j = 1:size(stats,2)
            corr(i,j) = stats{i,j}.corr;
            r2(i,j) = stats{i,j}.r2;
            mse(i,j) = stats{i,j}.mse;
            nmse(i,j) = stats{i,j}.nmse;
        end
    end
    
    p.corr = friedman(corr, 1, 'off');
    p.r2 = friedman(r2, 1, 'off');
    p.mse = friedman(mse, 1, 'off');
    p.nmse = friedman(nmse, 1, 'off');
    
    
else
    error('Unknown model type.')
end




end
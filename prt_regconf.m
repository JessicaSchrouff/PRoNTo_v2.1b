function out = prt_regconf(PRT, in,trainonly,d)
% function to remove confounds from raw data within 'in'
% trainonly is a flag set to 1 if removing using training data only
% otherwise do it using all data
% Copyright (C) 2016 Machine Learning & Neuroimaging Laboratory

% Written by A Rao
% $Id$

% copy input fields to output
out = in;

% get data and covs
traindata=in.train{d};
traincovs=in.tr_cov;
augtraincovs=[ones(size(traincovs,1)) traincovs];

if (trainonly==1)
    % adjust using training only
    out.train=traindata-augtraincovs*(pinv(augtraincovs)*traindata);
    %     out.test=testdata-augtestcovs*(pinv(augtraincovs)*traindata);
else
    
    testdata=in.test{d};
    testcovs=in.te_cov;
    augtestcovs=[ones(size(testcovs,1),1) testcovs];
    
    % adjust using all data
    alldata=[traindata;testdata];
    out.train=traindata-augtraincovs*(pinv([augtraincovs;augtestcovs])*alldata);
    out.test=testdata-augtestcovs*(pinv([augtraincovs;augtestcovs])*alldata);
end

end


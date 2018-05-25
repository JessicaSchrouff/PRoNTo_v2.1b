function [Ktr_r,Kte_r,Ktrte_r] = prt_remove_confounds_TrKernel(Ktr,Kte,Ktrte,Ctr,Cte)

% [Ktr_r, Ktr_te,Ktrte_r,R] = prt_remove_confounds_AR(Ktr,Kte,KtrKte,Ctr,Cte)
%
% Function to remove confounds from kernel.
% Ktr =Xtr*Xtr',Kte=Xte*Xte',Ktrte=Xtr*Xte'
% Ctr are the confounds (covariates) for training data
% Cte are the confounds for test data
% We return the deconfounded kernels as _r

%__________________________________________________________________________
% Copyright (C) 2015 Machine Learning & Neuroimaging Laboratory

% Written by A. Rao
% $Id$

% we need to include an intercept term to the confounds corresponding to a
% mean
% dont think its required if all kernels and confounds are centred
% if they are centred, it wont cause any harm having it (the mean should
% then be estimated as zero)

ntr=size(Ktr,1);
nte=size(Kte,1);

Ctr=[ones(ntr,1) Ctr]; % include intercept in confounds, can delete this line and one after if we tell users to include a column of ones in C_tr,C_te
Cte=[ones(nte,1) Cte];

% deconfound training kernel
Ctr_pi=pinv(Ctr); % pseudoinverse of training confounds
R=eye(ntr)-Ctr*Ctr_pi;
Ktr_r=R'*Ktr*R; % deconfounded training kernel

% deconfound test kernel
Ctetr_pi=Cte*Ctr_pi;
crossterm=Ctetr_pi*Ktrte;
Kte_r=Kte-crossterm-crossterm'+Ctetr_pi*Ktr*Ctetr_pi'; % deconfounded training kernel

% deconfound train test kernel
%Ktrte_r=-R*Ktr*Ctetr_plus'+R*Ktrte;
Ktrte_r=R*(-Ktr*Ctetr_pi'+Ktrte);
return
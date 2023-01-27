function pare = est_cal_like_hess(z,s,m,n,c,thi,restrict);
% est_cal_like_hess performs minimization of the quasi likelihood taking
% illconditioning of the Hessian into account. 
% Starting values are contained in thi.
%
% SYNTAX: pare = est_cal_like_hess(z,s,m,n,c,thi,restrict);
%
% INPUTS:  z ... Tx(s+m); observations
%          s ... integer; dimension of endogenous variables.
%          m ... integer; dimension of exogenous variables.
%          n ... integer; state dimension.
%          c ... integer; number of common trends.
%          thi ... theta structure containing the initial estimates.
%          restrict ... restrict structure used for calculating the system
%                  for given parameters. 
% 
% OUTPUTS:  pare ... dx1 optimal parameter value
%
% REMARKS: + optimization is performed repeatedly (5 iterations hard
% coded!) where in each iteration the parameters are reweighted by inverses
% of the square roots of the diagonal values of the estimated Hessian. This
% makes the problem more well conditioned. 
%         + values of rescaling are passed on to the parameterization in
%         restrict.scale. 
% 
% AUTHOR: dbauer, 14.1.2020.
parami = syst2param(thi,c,restrict);

% improve estimate
%options = optimoptions('fminunc','display','iter');
options = optimoptions('fminunc','display','iter','MaxIterations',1000);
options.MaxFunctionEvaluations = 20000;
%restrict.scale = ones(length(parami),1);
%restrict.Omega = thi.Omega; 

sizOm = s*(s+1)/2;
Pbull = -1; 

if isfield(restrict,'scale')
    restrict.scale = restrict.scale(sizOm+1:end);
end;

[pare,fval,exitflag] = fminunc(@(x) cal_quasi_like(x,z,s,m,n,c,Pbull,restrict),parami(sizOm+1:end),options);

for j=1:5
    restrict.scale = ones(length(pare),1);
    [llc2,resc2] = cal_quasi_like(pare,z,s,m,n,c,Pbull,restrict);
    Omegai = resc2'*resc2/size(resc2,1);
    %restrict.Omega = Omegai;
    sh = ones(length(pare),1);
    if length(pare)<20
        hess = hessian(@(x) cal_quasi_like(x,z,s,m,n,c,Pbull,restrict),pare);
        dh = max(0.001,diag(hess)); % regularize for negative values.
        sh = sqrt(dh);
    end
    if Pbull<0 
        restrict.scale = [1./sh];
    else
        restrict.scale = [ones(sizOm,1);1./sh];
    end
    pare = pare.*sh;
    [pare,fval,exitflag] = fminunc(@(x) cal_quasi_like(x,z,s,m,n,c,Pbull,restrict),pare,options);  
    pare = pare./sh;
end;

restrict.scale = ones(length(pare),1);
[llc2,resc2] = cal_quasi_like(pare,z,s,m,n,c,Pbull,restrict);
Omegai = resc2'*resc2/size(resc2,1);
    
% rescale
paraomi = extr_lowtri(Omegai);
pare = [paraomi;pare];



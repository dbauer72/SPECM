function result = compile_results(parc,s,m,n,c,restrict,y,dt,Pbull,nmax);
% compile_results writes all information corresponding to the parameter
% vector parc into the structure result. 
% 
% SYNTAX: result = compile_results(parc,s,m,n,c,restrict,y,dt,Pbull,nmax);
% 
% INPUT: parc ... parameter vector; is passed tp param2syst.
%        s    ... integer; dimension of endogenous vars
%        m    ... integer; dimension of exogenous vars
%        n    ... integer; state dimension
%        c    ... integer; number of common trends
%        restrict ... structure; contains information for parameterization.
%        y    ... Txs; observations of endogenous vars.
%        dt   ... Txm; observations of exogenous vars.
%        Pbull ... indicator; if Pbull>0: state started with stationary
%                  dist; else with zero.
%       nmax  ... maximal state dimension allowed (contained for
%       information in result.call. 
%
% OUTPUT: result: est_result structure (see there for description of
% fields.
%
% REMARKS: + system is calculated using param2syst; 
%          + quasi likelihood calculated using cal_quasi_like. 
%       
% AUTHOR: dbauer, 14.1.2020.
if isfield(restrict,'scale')
    restrict = rmfield(restrict,'scale');
end;
[~,~,~,~,~,thc] = param2syst(parc,s,n,m,c,restrict);
thc.B = zeros(n,m);

if Pbull<0 
    sizOm = s*(s+1)/2;
else
    sizOm = 0;
end;
[llc,resc] = cal_quasi_like(parc(sizOm+1:end),[y,dt],s,m,n,c,Pbull,restrict);
Omega= resc'*resc/size(resc,1);
th.Omega = Omega;

% calculate the variance matrices 
[V_theta,V_syst,~] = cal_quasi_like_V(parc,[y,dt],s,m,n,c,restrict);

% fill in everything into result structure
result = est_result();
result.urs = c;
if c>0
    result.ur = 'I(1)';
else
    result.ur = 'I(0)';
end

result.n = n;
result.s = s;
result.y = y;
result.dt = dt;
result.theta = thc;
result.deviance = llc;
result.aic = llc + 2*length(parc);
result.bic = llc + log(size(resc,1))*length(parc);
result.res = resc;
result.restrict = restrict;
result.Pbull = Pbull;
result.call = sprintf('SPECM(z,%d,%d,%d,%d,%d)',s,n,c,nmax,Pbull);
result.param = parc;
result.V_theta = V_theta;
result.V_syst = V_syst;

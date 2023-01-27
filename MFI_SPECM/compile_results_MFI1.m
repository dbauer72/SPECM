function result = compile_results_MFI1(parc,s,S,n,urs,restrict,y,dt,Pbull,nmax);
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
m = size(dt,2); 
[~,~,~,~,~,thc] = param2syst_MFI1(parc,s,n,urs,restrict); 
thc.B = zeros(n,m);

if Pbull<0 
    sizOm = s*(s+1)/2;
else
    sizOm = 0;
end;
[llc,resc] = cal_quasi_like_MFI1(parc,[y,dt],s,S,n,urs,Pbull,restrict);
Omega= resc'*resc/size(resc,1);
thc.Omega = Omega;
% fill in everything into result structure
result = est_result();
result.ur = 'MFI(1)';
result.urs = urs;
result.S = S;
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
result.call = sprintf('SPECM_MFI1(z,%d,%d,%d,%d,%d)',s,S,n,nmax,Pbull);
result.param = parc;
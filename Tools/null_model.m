function result = null_model(y,dt);
% provides the result structure for the null model of a simple 
% regression y_t = D u_t + e_t under assumption of iid for e_t.
%
% SYNTAX: result = null_model(y,dt);
%
% INPUT: y  ... Txs real matrix of observations
%        dt ... Txm real matrix of deterministics.
%
% OUTPUT: result ... result structure. 
%
% AUTHOR: dbauer, 28.8.2020

[T,s] = size(y); 
[~,m]= size(dt);
Dhat = (dt\y)';

resc = y - dt*Dhat';
Omega = resc'*resc/T; 


result = est_result();
result.urs = 0;
result.ur = 'I(0)';
result.n = 0;
result.s = s;
result.y = y;
result.dt = dt;

thc = theta_urs();
thc.which = 'SS';
thc.A = zeros(0,0);
thc.B = zeros(0,m);
thc.K = zeros(0,s);
thc.D= Dhat;
thc.C = zeros(s,0);
thc.Omega = Omega; 
thc.m=m;
thc.ur = '';
result.theta = thc;

% deviance 
llc = T*log(det(Omega))+s*T;
result.deviance = llc;
result.aic = llc + 2*s*m;
result.bic = llc + log(size(resc,1))*s*m;
result.res = resc;
rest.det_res = 0;
result.restrict = rest;
result.Pbull = 0;
result.call = ['null model'];

paraomi = extr_lowtri(Omega);

result.param = paraomi;
function [A,K,C,D,Omega,th] = param2syst_I2(param,s,n,m,ds,restrict); 
% takes the parameter vector, the integers and the restrictions and
% converts it to system plus innovation variance.
%
% SYNTAX: [A,K,C,D,Omega,th] = param2syst_I2(param,s,n,m,c,restrict); 
%
% INPUT: param ... dx1 parameter vector.
%        s     ... integer; dimension of endogenous vars.
%        n     ... integer; state dimension.
%        m     ... integer; dimension of exogenous vars.
%        ds     ... pair of integer; number of common trends.
%        restrict ... structure describing the restrictions. 
%
% OUTPUT: A,K,C,D ... state space system.
%         Omega   ... sxs innovation variance.
%         th      ... thetat structure corresponding to A,K,C,D.
%         gr_syst ... if wanted in the output contains the derivative of
%                       the system matrices with respect to the parameters. 
%                      Array: gr_syst(1).A, .K, .C, .D, .Omega contains the
%                      derivative of A,K,C,D, Omega respectively w.r.t. the
%                      first parameter, gr_syst(2).A, ... w.r.t. to the
%                      second param.
%
% REMARK: uses the parameterization corr. to the param paper. 
%
% AUTHOR: dbauer, 21.11.2019

det_res = restrict.det_res;
if isfield(restrict,'scale')
    scale = diag(restrict.scale); 
else
    scale = eye(length(param));
end;

param = scale*param;

% extract parameters for Omega
parom = param(1:(s*(s+1)/2));
nom = length(parom);
param(1:(s*(s+1)/2))=[];


Omega = fill_lowtri(parom,s);


% convert params to matrices
nd = s*m;
if det_res 
    nd = nd- (s-c);
end;


[th,param] = param2th_I2(param,s,n,ds,restrict);

A = th.A;
K = th.K;
C = th.C;

% parameters for deterministics. 
D = reshape(param,s,m);

th.A =A;
th.K =K;
th.C =C;
th.D =D;
th.Omega =Omega;
th.ur = 'I(2)';


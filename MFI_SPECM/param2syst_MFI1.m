function [A,K,C,D,Omega,th,gr_syst] = param2syst_MFI1(param,s,n,urs,restrict,Omega); 
% takes the parameter vector, the integers and the restrictions and
% converts it to system plus innovation variance.
%
% SYNTAX: [A,K,C,D,Omega,th,gr_syst] = param2syst_MFI1(param,s,S,n,urs,restrict,Omega); 
%
% INPUT: param ... dx1 parameter vector.
%        s     ... integer; dimension of endogenous vars.
%        n     ... integer; state dimension.
%        urs   ... unit root structure.
%        restrict ... structure describing the restrictions. 
%        Omega ... sxs innovation covariance matrix. 
%
% OUTPUT: A,K,C,D ... state space system.
%         Omega   ... sxs innovation covariance.
%         th      ... theta_urs structure corresponding to A,K,C,D.
%         gr_syst ... if wanted in the output contains the derivative of
%                       the system matrices with respect to the parameters. 
%                      Array: gr_syst(1).A, .K, .C, .D, .Omega contains the
%                      derivative of A,K,C,D, Omega respectively w.r.t. the
%                      first parameter, gr_syst(2).A, ... w.r.t. to the
%                      second param.
%
% REMARK: uses the parameterization corr. to the param paper. 
%
% AUTHOR: dbauer, 14.10.2020

det_res = restrict.det_res;
if isfield(restrict,'scale')
    scale = diag(restrict.scale); 
else
    scale = eye(length(param));
end;

param = scale*param;

% extract parameters for Omega
if nargin<6 % Omega not contained
    parom = param(1:(s*(s+1)/2));
    nom = length(parom);
    %param(1:(s*(s+1)/2))=[];
else
    parom = extr_lowtri(Omega);
end

if nargout>4
    np = length(param)+length(parom);
    [Omega,dOmega] = fill_lowtri(parom,s);
else
    Omega = fill_lowtri(parom,s);
end

% convert params to matrices

if nargout>6 % gradients wanted?
    nth = length(param);
    [th,param,dth] = param2th_MFI1(param,s,n,urs,restrict);
    nth = nth-length(param);    
else
    [th,param] = param2th_MFI1(param,s,n,urs,restrict);
end

m = floor(length(param)/s);
D = reshape(param(1:s*m),s,m);
th.D = D;
th.urs = urs;
th.ur = 'MFI(1)';
th.Omega = Omega;

A = real(th.A);
K = real(th.K);
C = real(th.C);



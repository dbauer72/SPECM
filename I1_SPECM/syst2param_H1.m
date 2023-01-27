function [par,restrict] = syst2param_H1(th,c,restrict);
% syst2param collects the parameters corresponding to the system subject to
% the restrictions formulated in restrict.
%
% SYNTAX: par = syst2param_H1(th,c,restrict);
%
% INPUT: th ... theta structure
%        c  ... number of common trends
%        restrict ... hypotheses 
% 
% OUTPUT: par ... vector of parameters.
%
% AUTHOR: dbauer, 21.11.2019. 

paraomi = extr_lowtri(th.Omega);
% parameters not taking restriction into account. 
parami = th2param(th,c,0);
% find parameters due to C1.
[s,n]=size(th.C);
nparc = s*c-c*(c+1)/2;
parami(1:nparc)=[];

% find parameters corresponding to H0. 
C1 = proj_H1(th.C(:,1:c),restrict);

par = ortho2par_H1(C1,restrict.H.H,restrict.H.t);
restrict.nparc = length(par);

if restrict.det_res % restrict the highest order deterministic term.    
    C1 = th.C(:,1:c);
    D = th.D; 
    pard = C1'*D(:,1);
    D2 = D(:,2:end);
    if size(D,2)>1
        pard = [pard(:);D2(:)];
    end;
    par = [paraomi(:);par(:);parami(:);pard];
else
    par = [paraomi(:);par(:);parami(:);th.D(:)];
end


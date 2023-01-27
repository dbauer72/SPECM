function par = syst2param(th,c,restrict);
% syst2param collects the parameters corresponding to the system subject to
% the restrictions formulated in restrict.
%
% SYNTAX: par = syst2param(th,c,restrict);
%
% INPUT: th ... theta structure
%        c  ... number of common trends
%        restrict ... hypotheses 
% 
% OUTPUT: par ... vector of parameters.
%
% AUTHOR: dbauer, 21.11.2019. 

paraomi = extr_lowtri(th.Omega);
if isfield(restrict,'syst2par')
    parami = feval(restrict.syst2par,th,c,restrict);
else 
    parami = th2param(th,c,1);
end;

if (restrict.det_res)&&(c>0) % restrict the highest order deterministic term.    
    C1 = th.C(:,1:c);
    D = th.D; 
    pard = C1'*D(:,1);
    D2 = D(:,2:end);
    if size(D,2)>1
        pard = [pard(:);D2(:)];
    end;
    par = [paraomi(:);parami(:);pard];
else
    par = [paraomi(:);parami(:);th.D(:)];
end
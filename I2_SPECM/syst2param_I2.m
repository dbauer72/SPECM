function par = syst2param_I2(th,ds,restrict);
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
parami = th2param_I2(th,ds,0);
par = [paraomi(:);parami(:);th.D(:)];

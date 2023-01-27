function [alpha_m,Va_m] = standardize_alpha(alpha,baralpha,Va);
% standardizes alpha by transforming it such that baralpha'*alpha=I_r.
%
% SYNTAX: [alpha_m,Va] = standardize_alpha(alpha,baralpha,Va);
%
% INPUT: alpha    ... sxr real matrix.
%        baralpha ... sxr real standardization matrix.
%        Va       ... sr x sr variance matrix.
%
% OUTPUT: alpha_m ... sxr standardized real matrix.
%         Va_m    ... sr x sr real transformed variance matrix.
%
% REMARK: standardization works as alpha*inv(baralpha'*alpha).
%
% AUTHOR: dbauer, 17.12.2021

[s,r]=size(alpha); 

iba = inv(baralpha'*alpha);
alpha_m = alpha*iba; 

% vectorize to calculate variance.
Tr = kron(eye(s) - alpha*iba'*baralpha',iba');
Va_m = Tr*Va*Tr';


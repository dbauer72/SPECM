function paro = extr_lowtri(Omega);
% extracts parameters as lower triangular factor of square matrix Omega.
%
% SYNTAX: paro = extr_lowtri(Omega);
%
% INPUTS: Omega sxs innovation covariance matrix.
%
% OUTPUTS: paro ... s(s+1)/2 x 1 parameter vector.
%
% REMARK: CC'=Omega for lower triangular matrix C is calculated. 
%   paro = vech(C).
%
% AUTHOR: dbauer, 29.10.2019

s =size(Omega,1);
C = chol(Omega,'lower'); % lower triangular Cholesky factor

paro =[];
for i=1:s
    paro = [paro(:);C(i:end,i)];
end;

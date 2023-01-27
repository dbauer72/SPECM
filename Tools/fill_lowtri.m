function [Omega,dOmega] = fill_lowtri(parom,s);
% fills in elements in lower triangular Cholesky factor and squares it.
%
% SYNTAX: [Omega,dOmega] = fill_lowtri(parom,s);
%
% INPUTS:   parom ... dx1 parameter vector 
%           s     ... integer, size of Omega.
%
% OUTPUTS:  Omega ... sxs innovation covariance matrix.
%           dOmega ... sxsxd array of derivatives w.r.t. entries of parom.
%           (calculated only, if second output argument included).
%
% REMARK: parom = vech(C) where Omega = CC'.
%
% AUTHOR: dbauer, 29.10.2019

LOm = zeros(s,s);
if nargout>1
    dOmega = zeros(s,s,length(parom));
    dLOm = zeros(s,s,length(parom));
    cur = 0;
end

dpar = eye(length(parom));

for i=1:s
    LOm(i:end,i)=parom(1:(s-i+1));
    if nargout>1
        dLOm(i:end,i,cur+[1:(s-i+1)])=eye(s-i+1);
        cur = cur + s-i+1;
    end
    parom(1:(s-i+1))=[];
end

Omega = LOm*LOm';
if nargout>1
    for j=1:size(dOmega,3)
        dOmega(:,:,j)= LOm*squeeze(dLOm(:,:,j))'+ squeeze(dLOm(:,:,j))*LOm';
    end
end

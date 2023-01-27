function param1 = mat2param(C1,K1);
% converts C1 and K1 matrices into parameters assuming in canonical form
%
% SYNTAX: param1 = mat2param(C1,K1);
%
% INPUT: C1 ... sxc matrix;
%        K1 ... cxs matrix;
%
% OUTPUT: param ... dx1 parameter vector.
%
% REMARK: assumes that C1 and K1 correspond to the generic neighborhood of
% the canonical form such that C1*C1 = I_c and K1 is p.u.t.
%
% AUTHOR: dbauer, 27.10.2019.

pc = ortho2par_LR(C1);

pk = [];
c = size(K1,1);

for j=1:c 
    pk = [pk,K1(j,j:end)];
end

param1 = [pc(:);pk(:)];
function De =  fill_D(param,C1,s,m,restrict);
% fills in parameters into D taking restrictions in restrict into account. 
%
% SYNTAX: De =  fill_D(param,C1,s,m,restrict);
%
% INPUT:  param ... parameter vector.
%         C1    ... matrix C1'C1 = I_c of loadings of common trends.
%         s     ... integer; dimension of endogenous vars.
%         m     ... integer; dimension of exogenous vars.
%         restrict ... structure indicating the restrictions.
%
% OUTPUT: De  ... sxm matrix D.
%
% REMARK: if restrict.det>0 -> restrict first column of De to lie in the
% column space of C1. 
% 
% AUTHOR: dbauer, 14.1.2020. 


c = size(C1,2);
if c>0
    [Q,R]= qr(C1(:,1:c));
else
    Q = eye(s);
end;

if (restrict.det_res)&&(c>0)
    De = zeros(s,m);
    pc = param(1:c);
    De(:,1)=C1(:,1:c)*pc(:);
    param(1:c)=[];
    if m>1
        De(:,2:end)= reshape(param,s,m-1);
    end
else
    De = reshape(param,s,m);
end
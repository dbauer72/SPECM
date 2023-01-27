function Ck = approx_equalC(C,urs);
% approx_equalC approximates an optimal common C.
%
% SYNTAX: Ck = approx_equalC(C,urs);
%
% INPUT: C ... matrix in real canonical form
%        urs .. unit root structure
%
% OUTPUT: Ck ... real sxck matrix of approximated factors. 
%
% AUTHOR: dbauer, 15.12.2021. 

[s,n]=size(C);
VCk = zeros(s,s);

cur = 0;
for j=1:size(urs,1)
    ck = urs(j,3);
    if urs(j,2)==0 % real unit root 
        Ck = C(:,cur+[1:ck]);
        VCk = VCk+ Ck * Ck';
        cur = cur+ck;
    else % complex root 
        Ck = C(:,cur+[1:ck])+sqrt(-1)*C(:,cur+ck+[1:ck]);
        VCk = VCk+ real(Ck * Ck');
        cur = cur+2*ck;
    end
end

% compute svd decomposition. 
[u,s,v]=svd(VCk); 

Ck = u(:,1:ck);


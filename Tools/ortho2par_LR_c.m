function [par,QL] = ortho2par_LR_c(C);
% ortho2par_LR extracts parameters for C'C=I_c according to the LR
% canonical form for complex matrix C.
%
% SYNTAX: par = ortho2par_LR_c(C);
%
% INPUT:  C ... sxc complex valued matrix, C0'C0=I_c.
%
% OUTPUT: par ... dx1 parameter vector.
%
% REMARK: delivers the parameters for a matrix C'C=I using the parameter values par
% such that C = Q_L[I;0]Q_R where Q_L and Q_R are products of Givens
% rotations.
%
% AUTHOR: dbauer, 14.10.2020%
[n,p] = size(C);
par = [];

% first transform from the right to lower diagonal in the heading pxp
% matrix.

% theta_R
for i=1:(p-1)
    for j=(i+1):p
        [pa,Qi] = pars_c(C(i,i),(C(i,j)));
        par=[par(:);pa(:)];
        C(:,[i,j])=C(:,[i,j])*conj(Qi');
    end
end

% now C is lower triangular.
% now start from the back.
QL = eye(n);

for i=p:-1:1
    for j=p:(n-1)
        [pa,Qi] = pars_c(C(i,i),C(j+1,i));
        par=[par(:);pa(:)];
        
        C([i,j+1],:)=Qi*C([i,j+1],:);
        QL([i,j+1],:)=Qi*QL([i,j+1],:);
    end;
end;

% finish
% diagonal matrix contains complex values.
cc = diag(C(1:p,:));
par = [par(:);angle(cc)];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% subroutine to extract parameters   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pa,Qi]= pars_c(C1,C2);

a = abs(C1); b = abs(C2);
if (a>0)
    phia = angle(C1);
else
    phia=0;
end;

if (b>0)
    phib = angle(C2);
else
    phib = 0;
end;

if (a>0)
    pa(1)=atan(b/a);
else
    if (b>0)
        pa(1)=pi/2;
    else
        pa(1)=0;
    end;
end;

pa(2)=phib-phia;
% generate Givens rotation:
i = sqrt(-1);
Qi = [cos(pa(1)),sin(pa(1))*exp(-i*pa(2));-sin(pa(1))*exp(i*pa(2)),cos(pa(1))];


end
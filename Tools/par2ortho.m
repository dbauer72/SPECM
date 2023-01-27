function C = par2ortho(par,n,p);
% par2ortho transforms the parameters in par into the orthonormal column C.
% Inverse to ortho2par.
%
% SYNTAX: C = par2ortho(par,n,p);
%
% INPUT:  par ... dx1 parameter vector.
%         n,p ... dimensions of C. 
%      
% OUTPUT: C ... sxc matrix, C'C=I_p.
%
% REMARK: delivers the parameters for a matrix C'C=I using the parameter values par
% such that C = Q[I;0] where Q is the product of Givens
% rotations. 
%
% AUTHOR: dbauer, 28.10.2003

C = [eye(p);zeros(n-p,p)];
for i=p:-1:1
    for j=(n-1):-1:i
        pc = par(end);
        par = par(1:end-1);
        Q= [cos(pc),-sin(pc);sin(pc),cos(pc)];
        C([i,j+1],:)=Q*C([i,j+1],:);
    end;
end;
%C = C(:,[p+1:end]);


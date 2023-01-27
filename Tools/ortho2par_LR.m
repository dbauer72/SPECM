function [par,QL] = ortho2par_LR(C);
% ortho2par_LR extracts parameters for C'C=I_c according to the LR
% canonical form. 
%
% SYNTAX: par = ortho2par_LR(C);
%
% INPUT:  C ... sxc matrix, C0'C0=I_c, C=H*C0.
%      
% OUTPUT: par ... dx1 parameter vector.
%
% REMARK: delivers the parameters for a matrix C'C=I using the parameter values par
% such that C = Q_L[I;0]Q_R where Q_L and Q_R are products of Givens
% rotations. 
%
% AUTHOR: dbauer, 20.11.2019%
[n,p] = size(C);
par = [];

% first transform from the right to lower diagonal in the heading pxp
% matrix. 

% theta_R
for i=1:(p-1)
    for j=(i+1):p
        pa=0;
         if C(i,i)>0
            pa =atan(C(i,j)/C(i,i));
            par(end+1)=pa;
  %          C([i,j+1],:)=[cos(pa),sin(pa);-sin(pa),cos(pa)]*C([i,j+1],:);
        elseif (C(i,i)==0)
            if(C(i,j)>0)
                pa = pi/2;
                par(end+1)=pa;
            elseif (C(i,j)<0)
                pa = -pi/2;
                par(end+1)=-pi/2;
            else
                pa=0;
                par(end+1)=0;
            end;
        elseif(C(i,i)<0)
            pa =atan(C(i,j)/C(i,i))+pi;
            par(end+1)=pa;
        end;
        C(:,[i,j])=C(:,[i,j])*[cos(pa),sin(pa);-sin(pa),cos(pa)]';
    end
end

% now C is lower triangular. 
% now start from the back. 
QL = eye(n); 

for i=p:-1:1
    for j=p:(n-1)
        if C(i,i)>0
            pa =atan(C(j+1,i)/C(i,i));
            par(end+1)=pa;
  %          C([i,j+1],:)=[cos(pa),sin(pa);-sin(pa),cos(pa)]*C([i,j+1],:);
        elseif (C(i,i)==0)
            if(C(j+1,i)>0)
                pa = pi/2;
                par(end+1)=pa;
            elseif (C(j+1,i)<0)
                pa = -pi/2;
                par(end+1)=-pi/2;
            else
                pa=0;
                par(end+1)=0;
            end;
        elseif(C(i,i)<0)
            pa =atan(C(j+1,i)/C(i,i))+pi;
            par(end+1)=pa;
        end;
        C([i,j+1],:)=[cos(pa),sin(pa);-sin(pa),cos(pa)]*C([i,j+1],:);
        QL([i,j+1],:)=[cos(pa),sin(pa);-sin(pa),cos(pa)]*QL([i,j+1],:);   
    end;
end;

% finish        
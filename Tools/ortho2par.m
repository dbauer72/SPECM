function par = ortho2par(C);
% ortho2par extracts parameters from an orthonormal column.
%
% SYNTAX: par = ortho2par(C);
%
% INPUT:  C ... sxc matrix, C'C=I_c.
%      
% OUTPUT: par ... dx1 parameter vector.
%
% REMARK: delivers the parameters for a matrix C'C=I using the parameter values par
% The matrix is parameterized using Givens rotations
%
% dbauer, 28.10.2003
[n,p] = size(C);
par = [];
for i=1:p
    for j=i:(n-1)
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
    end;
end;

        
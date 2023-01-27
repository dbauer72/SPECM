function [par,Ch] = ortho2par_H0(C,H,t);
% ortho2par_H0 extracts parameters for the restriction H_0: C = H*C0. 
%
% SYNTAX: par = ortho2par_H0(C,H,t);
%
% INPUT:  C ... sxc matrix, C0'C0=I_c, C=H*C0.
%         H ... sxt matrix, H'*H=I_t.
%         t ... integer, dimension of H. 
%      
% OUTPUT: par ... dx1 parameter vector.
%
% REMARK: delivers the parameters for a matrix C'C=I using the parameter values par
% The matrix is parameterized using Givens rotations
%
% Adaptation to ortho2par_H0: here the left and right params are separated as
% in the param paper. 
%
% implements restriction encoded in H. 
% AUTHOR: dbauer, 20.11.2019
[s,c] = size(C);
r = s-c; 

par = [];

% transform according to H. 
Ch = H'*C; 

% transform to diagonal matrix using right multiplication.

for i=1:c
    if i>(t-r)
        ic = i+r;
    else
        ic = i;
    end;
    for j=(i+1):c
        if Ch(ic,i)>0 
            pa =atan(Ch(ic,j)/Ch(ic,i));
            par(end+1)=pa;
        elseif (Ch(ic,i)==0)
            if (Ch(ic,j)>0)
                pa = pi/2;
                par(end+1)=pa;
            elseif (Ch(ic,j)<0)
                pa = -pi/2;
                par(end+1)=pa;  
            else
                pa = 0;
                par(end+1)=pa;
            end;
        elseif (Ch(ic,i)<0) 
            pa =atan(Ch(ic,j)/Ch(ic,i))+pi;
            par(end+1)=pa;
        end
        Ch(:,[i,j])=Ch(:,[i,j])*[cos(pa),sin(pa);-sin(pa),cos(pa)]';
    end
end

% now matrix should be standardized 
% extract parameters from the left 

% now start from the back. 
for i=(t-r):-1:1
    for j=r:-1:1
        if Ch(i,i)>0
            pa =atan(Ch(t-r+j,i)/Ch(i,i));
            par(end+1)=pa;
  %          C([i,j+1],:)=[cos(pa),sin(pa);-sin(pa),cos(pa)]*C([i,j+1],:);
        elseif (Ch(i,i)==0)
            if(Ch(j+r,i)>0)
                pa = pi/2;
                par(end+1)=pa;
            elseif (Ch(t-r+j,i)<0)
                pa = -pi/2;
                par(end+1)=-pi/2;
            else
                pa=0;
                par(end+1)=0;
            end;
        elseif(Ch(i,i)<0)
            pa =atan(Ch(t-r+j,i)/Ch(i,i))+pi;
            par(end+1)=pa;
        end;
        Ch([i,t-r+j],:)=[cos(pa),sin(pa);-sin(pa),cos(pa)]*Ch([i,t-r+j],:);
    end;
end;

% finish        
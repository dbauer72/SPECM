function [C,QL] = par2ortho_LR(par,n,p);
% par2ortho_LR transforms the parameters in par into the orthonormal column C.
% Inverse to ortho2par_LR. 
%
% SYNTAX: C = par2ortho_LR(par,s,c);
%
% INPUT:  par ... dx1 parameter vector.
%         s,c ... dimensions of C. 
%      
% OUTPUT: C ... sxc matrix, C'C=I_c. 
%
% REMARK: calculates C = Q_L [I_c;0]Q_R where Q_L=Q_L(par) and Q_R =
% Q_R(par). 

% AUTHOR: dbauer, 28.10.2003
C = [eye(p);zeros(n-p,p)];

% multiply using the entries of theta_L
QL = eye(n);
for i=1:p
    for j=(n-1):-1:p
        pc = par(end);
        par = par(1:end-1);
        Q= [cos(pc),-sin(pc);sin(pc),cos(pc)];
        C([i,j+1],:)=Q*C([i,j+1],:);
        QL([i,j+1],:)=Q*QL([i,j+1],:);        
    end;
end;


% multiply from the right using theta_R
for i=(p-1):-1:1
    for j=p:-1:(i+1)
        pc = par(end);
        par = par(1:end-1);
        Q= [cos(pc),-sin(pc);sin(pc),cos(pc)]';
        C(:,[i,j])=C(:,[i,j])*Q;
    end
end

%C = C(:,[p+1:end]);


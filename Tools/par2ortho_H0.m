function C = par2ortho_H0(par,s,c,H0);
% par2ortho_H0 transforms the parameters in par into the orthonormal column C.
% Inverse to ortho2par_H0. Implements the restriction H0: C = H*Ch;
%
% SYNTAX: C = par2ortho_H0(par,s,c,H0);
%
% INPUT:  par ... dx1 parameter vector.
%         s,c ... dimensions of C. 
%         H0  ... structure determining H0: H0.t... number of cols in H.
%                       H0.H: constraint such that C=  H*Ch. 
%      
% OUTPUT: C ... sxc matrix, C=  H*Ch. 
%
% AUTHOR: dbauer, 28.10.2003

t = H0.t;
H = H0.H;

r = s-c;
RL = [eye(t-r);zeros(r,t-r)];
for i=1:(t-r)
    for j=1:r
        % rotate entries.
        pc = par(end);
        par = par(1:end-1);
        Q= [cos(pc),-sin(pc);sin(pc),cos(pc)];
        RL([i,t-r+j],:)=Q*RL([i,t-r+j],:);
    end
end

Ch = [RL,zeros(t,s-t);zeros(s-t,t-r),eye(s-t)];

for i=c:-1:1
    for j=c:-1:(i+1)
        pc = par(end);
        par = par(1:end-1);
        Q= [cos(pc),-sin(pc);sin(pc),cos(pc)]';
        Ch(:,[i,j])=Ch(:,[i,j])*Q;
    end
end

C = H*Ch;
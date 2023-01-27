function C = par2ortho_H1(par,s,c,H0);
% par2ortho_H1 transforms the parameters in par into the orthonormal column C.
% Inverse to ortho2par_H1. Implements the restriction H0: b'*C = 0;
%
% SYNTAX: C = par2ortho_H1(par,s,c,H0);
%
% INPUT:  par ... dx1 parameter vector.
%         s,c ... dimensions of C. 
%         H0  ... structure determining H0: H0.t... number of cols in H.
%                       H0.H: constraint such that H(:,t+1:end)'*C=0. 
%      
% OUTPUT: C ... sxc matrix.
%
% AUTHOR: dbauer, 28.10.2003

H = H0.H;
t = H0.t;

Ch = par2ortho(par,s-t,c);
C = H(:,t+1:end)*Ch;
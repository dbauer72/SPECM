function [t,d] = seig(A);
% returns the eigenvalues sorted accordingto their modulus
%
% SYNTAX: [t,d] = seig(A);
%
% INPUT: A ... nxn real matrix;
%
% OUTPUT: t ... nxn real matrix of eigenvectors.
%         d ... nxn diagonal matrix of eigenvalues.
%
% REMARK: eigenvalues are sorted descending corresponding to their absolute
%         values. 
%         A = t*d*inv(t).
%
% AUTHOR: dbauer, 19.8.2020.
[t,d] = eig(A);
[x,i] = sort(abs(diag(d)),'descend');
d = d(i,:);d=d(:,i);
t = t(:,i);
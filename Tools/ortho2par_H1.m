function [par] = ortho2par_H1(C,H,t);
% ortho2par_H1 extracts parameters for the restriction H_1: H(:,1:t)'C =0. 
%
% SYNTAX: par = ortho2par_H1(C,H,t);
%
% INPUT:  C ... sxc matrix, Ch'Ch=I_c, H'*C=0.
%         H ... sxt matrix, H'*H=I_t.
%         t ... integer, dimension of H. 
%      
% OUTPUT: par ... dx1 parameter vector.
%
% REMARK: delivers the parameters for a matrix C'C=I using the parameter values par
% The matrix is parameterized using Givens rotations
%
% Adaptation to ortho2par_H1: H(:,1:t)'C =0.
% in the param paper. 
%
% implements restriction encoded in H. 
% AUTHOR: dbauer, 20.11.2019%

[s,c] = size(C);

par = [];

% transform according to H. 
Ch = H(:,t+1:end)'*C; 

par = ortho2par(Ch);
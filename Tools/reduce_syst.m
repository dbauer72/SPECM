function [Ab,Bb,Cb,red_norm_hank] = reduce_syst(A,B,C,nr);
% reduce_syst provides a reduced order system achieved by  
% SVD of the Hankel matrix and setting singular values to zero. 
%
% SYNTAX: [Ab,Bb,Cb,red_norm_hank] = reduce_syst(A,B,C,nr);
%
% INPUT:  (A,B,C) ... state space system
%         nr ... integer; reduced order. 
%
% OUTPUT: (Ab,Bb,Cb) ... reduced order system
%         red_norm_hank ... real; norm of difference in Hankel matrix. 
%
% AUTHOR: dbauer, 4.9.2020.

[n,s]= size(B); 

H = myhank(A,B,C,3*n);

[U,S,V]= svd(H);

Of = U(:,1:nr);
Cf = S(1:nr,1:nr)*V(:,1:nr)';
Bb = Cf(:,1:s);
Cb = Of(1:s,:);
Ab = Of(1:(end-s),:)\Of(s+1:end,:);

Hre = myhank(Ab,Bb,Cb,3*n);

red_norm_hank = norm(H-Hre);

function [th,Ae,Ke,Ce] = ss2ech_n(A,K,C);
% transfer the transfer function into (reordered) echelon canonical form
% such that On = [I_n;..].
%
% SYNTAX: [th,Ae,Ke,Ce] = ss2ech(A,K,C);
%
% INPUT: A,K,C ... triple of system matrices.
%       
% OUTPUT: 
%          th ... theta structure corresponding to A,K,C.
%          Ae,Ke,Ce ... triple of system matrices in reordered can. form.
% 
% COMMENT: assumes generic neighborhood.
%
% AUTHOR: dbauer, 27.10.2019

[n,s] = size(K);

On = zeros(n*s,n); 
Cn = zeros(n,n*s);

On(1:s,:)=C;

for j=2:n
    On((j-1)*s+[1:s],:)=On((j-2)*s+[1:s],:)*A;
end

% transformation then uses the first n rows of On.
Tr = On(1:n,:);
if rank(Tr)<n
    Tr = eye(n);
end

iTr= inv(Tr);
% transform 
Ae = Tr*A*iTr;
Ke = Tr*K;
Ce = C*iTr;

% convert to real, if necessary: 
Ae = real(Ae);
Ke = real(Ke);
Ce = real(Ce);

% fill into theta object.
th = theta_urs;
th.which = 'SS';
th.A = Ae;
th.C= Ce;
th.K = Ke;
th.Omega = eye(s);

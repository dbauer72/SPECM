function parambull = mat2param_bull(Abull,Kbull,Cbull);
% converts matrices in echelon canonical form into parameter vector.
%
% SYNTAX: parambull = mat2param_bull(Abull,Kbull,Cbull);
%
% INPUT: A,K,C ... triple of system matrices.
%
% OUTPUT: param ... dx1 vector of parameters. 
%
% REMARK: extracts the parameters for systems in particular canonical form:
%   For s>= n: C = [I_n;Co], A and K only contain free parameters. 
%  For s<n: C = [I_s,0], A = [0,I_s;Ao], K only contains free entries. 
%
% AUTHOR: dbauer, 27.10.2019.

[n,s]=size(Kbull);

% if s >= n then parameters are contained in C, A and K are full of
% parameters. 

if (s>= n)
    Ch = Cbull(n+1:end,:);
    parambull = [Ch(:);Abull(:);Kbull(:)];
else
    % things are more complicated here. 
    % no parameters in C,
    % some entries in A are constrained, Kbull fully parameterized. 
    Ah = Abull(n-s+1:end,:);
    parambull = [Ah(:);Kbull(:)];
end

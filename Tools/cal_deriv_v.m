function [dv] = cal_deriv_v(dy,Abar,K,tilC);
% cal_deriv_v derives the filtered state k_-(L)delta y_t w.r.to. the
% parameters contained in (Abar,K,C) in the generic neighborhood where O
% starts with the identity matrix and hence all parameters are in the last s rows of Abar
% and all of K. 
% 
% SYNTAX: [dv] = cal_deriv_v(dy,Abar,K,C);
%
% INUPT:  dy ... Txs real matrix of delta y_t. 
%         Abar ... nxn real matrix of inverse tf, Abar = A-KC
%         K   ... nxs real matrix.
%         tilC ... C(I-Abar)^{-1}Abar. 
%
% OUTPUT: dv ... cell vector of dimension s of Tx (2ns) matrix of
% derivatives.
%
% REMARK: uses the C(alpha) approach to cope with estimation of hat v_t.
%
% AUTHOR: dbauer, 25.10.2021.

[n,s] = size(K,1);
v0 = zeros(n,1); 

vt = ltitr(Abar,K,dy,v0);

dv = cell(s);
Is = eye(s); 

for j=1:(n*s)
    for cj = 1:n
        dAvt = vt(:,cj);
        for rj = 1:s 
            dvj = ltitr(Abar,Is(:,rj),vt);
        end
    end
    
end


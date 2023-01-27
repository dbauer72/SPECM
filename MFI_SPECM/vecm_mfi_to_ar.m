function ar = vecm_mfi_to_ar(est,S);
% rewrites the VECM MFI1 of Johansen Schaumburg into plain VAR in order to
% obtain predictions 
%
% SYNTAX: ar = vecm_mfi_to_ar(est);
%
% INPUT: est ... structure for MFI(1) estimation.
%        S   ... integer; number of observations per year. 
%
% OUTPUT: ar ... s x s*(k+S) matrix of AR coefficients. 
%
% REMARK: 
% format: y_t = a_1 y_{t-1} + ... + a_k y_{t-k} + e_t
% where ar = [a_1,a_2,...,a_k]. 
%
% AUTHOR: dbauer, 19.8.2020 


Gamma = est.Gamma;

[s,ps]=size(Gamma);
p = ps/s;

% reorder
%reord = [1:k;k+1:2*k];
%Gamma = Gamma(:,reord(:));


ar = zeros(s,(p+S)*s); % add four lags due to 1-L^S. 
% stationary coefficients due to Gamma. 
ar(:,1:ps)=Gamma; 
ar(:,s*S+[1:ps])=ar(:,s*S+[1:ps])-Gamma;
% add S-th lag.
ar(:,(S-1)*s+[1:s])=ar(:,(S-1)*s+[1:s])+eye(s);

M = length(est.beta); % number of unit roots possible. 
% cycle over roots 
for j=1:M
    if (j>1)&&(j<M) % complex root 
        fr = exp(2*pi*sqrt(-1)*(j-1)/S);
        polm = [1,conj(fr)];
        polm(end+1)=0;
        polm = polm - [0,polm(1:end-1)];
        polm(end+1)=0;
        polm = polm + [0,polm(1:end-1)];
        for k=2:M-1
            if (k ~= j) 
                frk =  exp(2*pi*sqrt(-1)*(k-1)/S);
                polm(end+1) = 0;
                polm = polm - frk*[0,polm(1:end-1)];
                polm(end+1) = 0;
                polm = polm - conj(frk)*[0,polm(1:end-1)];                
            end
        end
        ah = est.alpha{j}/2; 
        rah = size(ah,2)/2;
        alpha = ah(:,1:rah)-sqrt(-1)*ah(:,rah+[1:rah]);
        bh = est.beta{j};         
        beta = bh(1:s,1:rah)-sqrt(-1)*bh(1:s,rah+[1:rah]);
        
        Pi = alpha*beta';
        
        polm =polm/(polm*(conj(fr).^([0:S-1]))')/fr;
        
        for k=1:S 
            ar(:,(k-1)*s+[1:s])=ar(:,(k-1)*s+[1:s]) + 2*real(Pi(:,1:s)*polm(k)); 
        end
    end
    if j==1 % root z=1.
        pz = ones(1,S)/S;
        Pi = est.alpha{j}*est.beta{j}'; 
        for k=1:S 
            ar(:,(k-1)*s+[1:s])=ar(:,(k-1)*s+[1:s]) + Pi(:,1:s)*pz(k); 
        end
    end
    if j==M % roots z=-1 
        pz = (-1).^[1:S]/S;
        Pi = est.alpha{j}*est.beta{j}'; 
        for k=1:S 
            ar(:,(k-1)*s+[1:s])=ar(:,(k-1)*s+[1:s]) + Pi(:,1:s)*pz(k); 
        end
    end
end

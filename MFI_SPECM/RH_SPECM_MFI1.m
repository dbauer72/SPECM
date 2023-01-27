function [LL] = RH_SPECM_MFI1(y,S,Abar,B,Joh_j);
% Ribarits_hanzon idea for estimation of RH_VECM_MFI1 form in the MFI1
% case.
%
% SYNTAX: [LL] = RH_VECM(y,Abar,B,r,rest,Joh_j);
%
% INPUT: y    ... Txs data matrix,
%        S    ... integer; number of seasons;
%        Abar ... nxn matrix A-BC.
%        B    ... nxs matrix
%        Joh_j ... integer; indicating the version of Johansen's
%        deterministics. 
%
% OUTPUT: LL   ... minimized log likelihood.
%
% REMARK: the procedure implements the 'un-restricted' estimate regressing
% out the other unit root components. The likelihoods for each unit root
% are calculated independent of the other roots. 
%
% AUTHOR: dbauer, 22.10.2020. 

[T,s]=size(y); 
% prepare regressions 
if nargin< 5
    Joh_j = 5; % no deterministic terms. 
end;


n = size(Abar,1);
% steps: calculate tilde x_t: filtered version of Delta y.
[pX,pXk,Xm2,~,urs_full] = cal_reg_urs(y,S,s); 

% deal with deterministics 
[pXk,Xm2] = cal_reg_det(S,pXk,Xm2,urs_full,Joh_j);

pXk = zeros(size(pXk,1),0); % for the state space case, different stat. regressors are needed. 

% no rename to conform to the notation of Lukas Matuschek.
Z0t = pX;
v0 = zeros(n,1);
tB = ((-inv(eye(n)-Abar^S)*Abar^S)*B);
vt = ltitr(Abar,tB,pX,v0);

Z2t = vt;

% iterate over all roots
for j=1:size(urs_full,1)
    pXk = Z2t;
    for k=1:size(urs_full,1) % cycle over unit roots
        if (k ~= j)
            XX = squeeze(Xm2(:,:,k));
            if urs_full(k,2)==1 % complex unit root
                pXk = [pXk,XX];
            else % real unit root
                pXk = [pXk,XX(:,1:s)];
                if Joh_j == 4 % add regressor to coint. relation.
                    pXk = [pXk,XX(:,2*s+1)]; 
                end
                if (Joh_j<4)&&(urs_full(j,1)>0.0001)
                    pXk = [pXk,XX(:,2*s+1)]; 
                end
            end
        else
            XX = squeeze(Xm2(:,:,k));
            if urs_full(j,2)==1 % complex unit root
                Z1t = XX(:,1:2*s);
            else
                Z1t = XX(:,1:s);
            end
        end    
    end
    
    % now extract all regressors 
    R0t = Z0t - pXk*(pXk\Z0t);
    R1t = Z1t - pXk*(pXk\Z1t);
    
    if urs_full(j,2) == 1 % complex unit root
        Rit = R1t(:,1:s)+sqrt(-1)*R1t(:,s+[1:s]); % complex valued X_t^{(m)}
        if size(R1t,2)>2*s % determinstics contained?
            Rit = [Rit,R1t(:,2*s+1)+sqrt(-1)*R1t(:,2*s+2)];
        end;
        
        S00 = R0t'*R0t/(T-S);
        S11 = R1t'*R1t/(T-S);
        S01 = R0t'*R1t/(T-S);
        LR(1) = real((log(det(S00))));
        for r=1:s
            [~,LR(r+1)] = cal_betam(S00,S01,S11,r);
        end
        LL{j}=(-LR)*(T-S);
    else
        S00 = R0t'*R0t/(T-S);
        S11 = R1t'*R1t/(T-S);
        S01 = R0t'*R1t/(T-S);               
        S11d0 = S11 - S01'*(S00(1:s,1:s)\S01);
        S11d0(isnan(S11d0))=0;
        S11d0(isinf(S11d0))=0;
        me = min(eig(S11));
        if me<10^(-4)
            S11 = S11+eye(s)*(10^(-4)-me);
            S11d0 = S11d0+eye(s)*(10^(-4)-me);
        end;
        me = min(eig(S00));
        if me<10^(-4)
            S00 = S00+eye(s)*(10^(-4)-me);
        end;
        me = min(eig(S11d0));
        if me<10^(-4)
            S11d0 = S11d0+eye(s)*(10^(-4)-me);
        end;
                
        c11 = chol(S11);
        c00= chol(S00);
        B01 = inv(c00)'*S01*inv(c11);
        [~,~,V] = svd(B01);
        
        LR(1) = real(-T*(log(det(S00))));
        for r=1:s
            Betam = inv(c11)*V(:,1:r);
            LR(r+1) = real(-(T-S)*(log(det(S00))+log(det(Betam'*S11d0*Betam)/det(Betam'*S11*Betam))));
        end
        %LR = LR-LR(1);
        LL{j}= LR;
     end
end



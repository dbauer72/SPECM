function [LL,th_est,res,LR,estimate,u] = VECM_MFI1(y,p,S,urs,Joh_det);
% VECM_MFI1 uses the Johansen Schaumburg setting to estimate the VECM form
% for the MFI(1) model (seasonal cointegration),
% in a VAR Model in Error Correction Format for seasonal unit roots:
% p(L) X_t = sum a_mb_m' X_{t-1}^{(m)} + G1*p(L)X_{t-1} + ... + Gk*p(L)X_{t-k} + e_t
%
% SYNTAX:   [LL,th_est,res] = VECM_MFI1(y,p,S,urs,Joh_det);
%
% Input:     y ... Txs real; observations, 
%            p ... integer; number of Lagged difference terms in VECM
%                       formulation 
%            S ... integer; number of seasons.
%          urs ... Mx3 matrix. information on unit roots. j-th row:
%                       [omega_j,compl?,c_j].
% Output:    
%        LL      ... real; log of cond. quasi likelihood value
%        th_est  ... theta structure for estimated VAR system.
%
% This procedure utilizes the procedure to solve the
% generalized eigenvalue problem provided by J.P. LeSage,
% which has been modified to satisfy our needs, JPL_JOHA!
%
% NOTE: No deterministic variables included in this version!
% EXTERNAL FUNCTIONS: trimr, lag, tdiff, jpl_joha
%-------------------------------------------------------------
% function adapted from original procedures by M. Wagner. 
% dbauer, 21.12.2015; 19.8.2020

% entry checks
[T,s] = size(y);

if T<s
    error('observations need to be provided in format Txs!');
end;

[M,~] = size(urs);

if isempty(p)
    p = -fix(sqrt(T));
    disp('No lag length supplied, estimating using AIC up to sqrt(T)!');
end;

if p<0 
    p = aicest(y,s,-p);
end;

% estimation procedure 
% generate regressors 
% generates always all regressors, even if only fewer unit roots are
% contained. 
[pX,pXk,Xm2,~,urs_full] = cal_reg_urs(y,S,p); 

% deal with deterministics 
[pXk,Xm2] = cal_reg_det(S,pXk,Xm2,urs_full,Joh_det);

% find out, which unit roots are present. 
urs(:,4)=0;
for j=1:size(urs_full,1)
    abu = abs(urs(:,1)-urs_full(j,1));
    [mi,mind]=min(abu);
    if mi>0.0001 %unit root not contained
        XX = squeeze(Xm2(:,:,j));
        if urs_full(j,2)==1
            pXk = [pXk,XX];
        else
            pXk = [pXk,XX(:,1:s)];
            if Joh_det == 4 % add regressor to coint. relation.
                pXk = [pXk,XX(:,2*s+1)]; 
            end;
            if (Joh_det<4)&&(urs_full(j,1)>0.0001)
                pXk = [pXk,XX(:,2*s+1)]; 
            end;
        end
    else % unit root contained -> remember index. 
        urs(mind,4)=j;
    end
end


in = find(urs(:,4)>0); % only these roots are contained. 



% initial estimate 
estimate = initial_sc(pX,S,p,pXk,Xm2,urs,Joh_det);

% concentrate out pXk.
R0t = pX - pXk * (pXk\pX);
for j=1:length(in) 
    jj = urs(in(j),4);
     XX = squeeze(Xm2(:,:,jj));
     R1t(:,:,j) = XX - pXk * (pXk\XX);
end;

ok = 1;
% iterate over unit roots.
while ok<3
    for m=1:M % cycle over roots
        % estimate alpha and beta
        [estimate,LR{m}] = improve_alphabeta(R0t,R1t,estimate,m,'subtract',Joh_det);
        
        % update results
        [estimate,u,Dhat] = update_estimate(estimate,pX,pXk,Xm2,Joh_det);
    end;
    ok = ok+1;
end;

% now generate the theta structure. 
ar = vecm_mfi_to_ar(estimate,S);
th_est = theta_urs();

th_est.which = 'poly';
th_est.a = [eye(s),-ar];
th_est.b= eye(s);
th_est.d = Dhat;

th_est.ur = 'MFI(1)';
th_est.urs = urs;

% generate residuals.
Z = u(S+1:end,:);
for j=1:(p+S)
    Z = [Z,y(p+S-j+[1:(T-p-S)],:)];
end;


% calculate residuals
res = y(p+S+1:end,:) - Z*[Dhat,ar]';
Omega = res'*res/(T-p-S);
th_est.Omega = Omega; 

LL = -log(det(Omega))*(T-p-S)/2 - 0.5*(T-p-S)*s;

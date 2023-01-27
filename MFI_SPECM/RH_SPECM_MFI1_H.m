function [LL,df] = RH_SPECM_MFI1_H(y,S,Abar,B,Joh_j,urs,restrict)
% Ribarits_hanzon idea for estimation of RH_VECM_MFI1 form in the MFI1
% case.
%
% SYNTAX: [LL] = RH_SPECM_MFI1_H(y,Abar,B,r,rest,Joh_j,urs,restrict);
%
% INPUT: y    ... Txs data matrix,
%        S    ... integer; number of seasons;
%        Abar ... nxn matrix A-BC.
%        B    ... nxs matrix
%        Joh_j ... integer; indicating the version of Johansen's
%        deterministics. 
%        urs  ... unit root structure providing integers of common cycles.
%        restrict ... .type: 'A' real beta, 'B' real alpha, 'C' both. 
%                     .beta: frequencies where beta is restricted to be
%                     real.
%                     .alpha: frequencies where alpha is restricted to be
%                     real. 
%
% OUTPUT: LL   ... minimized log likelihood.
%
% REMARK: + The procedure implements the 'un-restricted' estimate regressing
% out the other unit root components. The likelihoods for each unit root
% are calculated independent of the other roots. 
%
%  + real beta is enforced in the maximization of cal_betam_R. 
%  + real alpha leads to estimation as in real unit root case: SVD of S_01
%  S_00^(-1) S_01 S_11^(-1).
%
% AUTHOR: dbauer, 15.2.2022. 

[T,s]=size(y); 
% prepare regressions 
if nargin< 5
    Joh_j = 5; % no deterministic terms. 
end;


n = size(Abar,1);
% steps: calculate tilde x_t: filtered version of Delta y.
[pX,pXk,Xm2,~,urs_full] = cal_reg_urs(y,S,s); 

% find the frequencies with restrictions. 
res_beta = zeros(size(urs_full,1),1);
res_alpha = zeros(size(urs_full,1),1);
df = res_beta;

switch restrict.type 
    case 'A'
        for j=1:length(restrict.betafr)
            inf = find(abs(urs_full(:,1)-restrict.beta(j))<0.001);
            if ~isempty(inf)
                res_beta(inf)=1;
            end;
        end
    case 'B' 
        for j=1:length(restrict.alphafr)
            inf = find(abs(urs_full(:,1)-restrict.alpha(j))<0.001);
            if ~isempty(inf)
                res_alpha(inf)=1;
            end
        end
    case 'C'
        for j=1:length(restrict.betafr)
            inf = find(abs(urs_full(:,1)-restrict.beta(j))<0.001);
            if ~isempty(inf)
                res_beta(inf)=1;
            end
        end
        for j=1:length(restrict.alphafr)
            inf = find(abs(urs_full(:,1)-restrict.alpha(j))<0.001);
            if ~isempty(inf)
                res_alpha(inf)=1;
            end
        end
end

urs_full(:,3)=0;
for j=1:size(urs,1)
    inf = find(abs(urs_full(:,1)-urs(j,1))<0.001);
    if ~isempty(inf)
        urs_full(inf,3)=urs(j,3);
    end
end

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
        if res_beta(j) ==0     % no restriction on beta       
            if (res_alpha(j)==0) % no restriction on alpha
                [~,LR] = cal_betam(S00,S01,S11,s-urs_full(j,3));
            else % alpha restricted to be real
                Betam = real_SVD(S00,S01,S11,s-urs_full(j,3));
                Rbeta = R1t*Betam;
                alpham = (Rbeta\R0t)';
                res_r = R0t- Rbeta*alpham';
                LR = log(det(res_r'*res_r/(T-S)));
                df(j) = (s-urs_full(j,3))*urs_full(j,3); % dfs: alpha' = [I,alpha_r] where alpha_r is real -> number of entries in alpha_r restricted.
            end
        else % beta restricted to be real
             if (res_alpha(j)==0) % alpha not restricted
                [~,LR] = cal_betam_R(S00,S01,S11,s-urs_full(j,3));
                df(j)= (s-urs_full(j,3))*urs_full(j,3); % dfs: beta' = [I,beta_r] where beta_r is real -> number of entries in beta_r restricted. 
             else % both are restricted.
                Betam = real_SVD(S00,S01(:,1:s),S11(1:s,1:s),s-urs_full(j,3));
                Rbeta = R1t(:,1:s)*Betam;
                alpham = (Rbeta\R0t)';
                res_r = R0t- Rbeta*alpham';
                LR = log(det(res_r'*res_r/(T-S)));
                df(j) = (2*s-urs_full(j,3))*urs_full(j,3); % one of alpha and beta normalized such that only the lower part is restricted. The other one is the fully restricted. 
            end
       end
        LL{j}=(-LR)*(T-S);
    else
        if (~isempty(strmatch(restrict.type(1:4),'I1.1')))&&(urs_full(j,1)==0) % linear restrictions at z=1.
            rest = restrict;
            rest.type = restrict.type(6:end);
            r = s-urs_full(j,3);
            fprintf('\n Testing for linear restrictions on alpha_0 beta_0^T at z=1 \n');
            [LR,alphahat,betahat,df(j)] = RH_specm_H0_Sij(R0t,R1t,r,Joh_j,'n',rest);
        elseif (~isempty(strmatch(restrict.type(1:5),'I1.m1')))&&(urs_full(j,1)==1) % linear restrictions at z=-1.
            rest = restrict;
            rest.type = restrict.type(7:end);
            r = s-urs_full(j,3);
            fprintf(' \n Testing for linear restrictions on alpha_0 beta_0^T at z=-1 \n');
            [LR,alphahat,betahat,df(j)] = RH_specm_H0_Sij(R0t,R1t,r,5,'n',rest);
        else

            S00 = R0t'*R0t/(T-S);
            S11 = R1t'*R1t/(T-S);
            S01 = R0t'*R1t/(T-S);               
                
            S11d0 = S11 - S01'*(S00(1:s,1:s)\S01);
            Betam = real_SVD(S00,S01,S11,s-urs_full(j,3));
            LR = real(-(T-S)*(log(det(S00))+log(det(Betam'*S11d0*Betam)/det(Betam'*S11*Betam))));
        end
        LL{j}= LR;
     end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  estimate Beta using real SVD calculations     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Betam = real_SVD(S00,S01,S11,r);

s = size(S00,1);

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

Betam = inv(c11)*V(:,1:r);


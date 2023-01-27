function [LL,alphahat,betahat,df] = VECM_H0(y,k,r,Joh_j,restrict)
% VECM formulation for testing restrictions on beta 
%
% SYNTAX: [LL,alphahat,betahat,df] = VECM_H0(y,k,r,Joh_j,restrict);
%
% INPUT: y    ... Txs data matrix,
%        k    ... lag order in the VECM,
%        r    ... rank of alpha beta' 
%        Joh_j ... integer; indicating the version of Johansen's
%        deterministics.
%        restrict ... specifies H0.
%
% OUTPUT: LL   ... 2x1 minimized log likelihood under H0 and (H0 u H1).
%         alphahat, betahat ... sxr matrices under the restriction
%         df  ... integer, degrees of freedom.
%
% COMMENT:  + implements estimation for three different restrictions:
%           (A) beta = restrict.beta.
%           (B) beta = restrict.H* phi
%           (C) beta = [restrict.b,beta_2]
%           restrict.type = "A"|"B"|"C".
%          + for Joh_j = 2,4 the additional regressor in beta is included
%          at the head not at the tail!
%
% AUTHOR: dbauer, 2.5.2022.

LL = [0,0];
% prepare regressions
if nargin< 4
    Joh_j = 5; % no deterministic terms.
end

if nargin <5
    error('VECM_H0: Restrictions need to be specified. See help for details.');
end

[T,s] = size(y);

Z0t = diff(y);
Z1t = y(1:end-1,:);
Z2t = zeros(T-1,s*k);

for j=1:k
    Z2t(:,(j-1)*s+[1:s])=[NaN(j,s);Z0t(1:end-j,:)];
end;

Z0t = Z0t(k+1:end,:);
Z1t = Z1t(k+1:end,:);
Z2t = Z2t(k+1:end,:);

Z2to= Z2t;
% here comes the separation corresponding to restriction or not. 
plus1= 0;
dt = zeros(T,0);

Z1to = Z1t;

switch Joh_j
    case 1
        Z0t = detrend(Z0t,1);
        Z1t = detrend(Z1t,1);
        Z2t = detrend(Z2t,1);
        dt = [ones(T,1),[1:T]'];
    case 2
        Z0t = detrend(Z0t,0);
        Z1t = [[1:T-k-1]',Z1t];
        Z1to = Z1t;
        Z1t = detrend(Z1t,0);
        Z2t = detrend(Z2t,0);
        dt = [ones(T,1)];
        plus1=1;
    case 3
        Z0t = detrend(Z0t,0);
        Z1t = detrend(Z1t,0);
        Z2t = detrend(Z2t,0);
        dt = [ones(T,1)];
    case 4
        Z1t = [ones(T-k-1,1),Z1t];
        Z1to = Z1t;
        plus1=1;
end

Phi0 = (Z2t\Z0t);
Phi1 = (Z2t\Z1t);

% generate R_{0t}
R0t = Z0t - Z2t*Phi0;

% generate R_{1t}
R1t = Z1t - Z2t*Phi1;

% estimation without taking restrictions into account
[alphahat,betahat,~,~]= est_alpha_beta(R0t,R1t,r,Joh_j);
% unrestricted residuals
res_u = R0t - R1t*betahat*alphahat';

LL(1) =  (T-2)*log(det(res_u'*res_u/(T-2)))+(T-2)*s;

% estimation under H0.
switch restrict.type
    case 'A' % betahat = beta.
        if ~isfield(restrict,'beta')
            error('VECM_H0: Type A hypotheses needs to specify beta!');
        end
        betahat = restrict.beta;
        %r = size(betahat,2);
        R1beta = R1t * betahat;
        alphahat = (R1beta\R0t)';
        res_r = R0t - R1beta*alphahat';
        LL(2) =  (T-2)*log(det(res_r'*res_r/(T-2)))+(T-2)*s;
        df = r*(s+plus1-r);
    case 'B' % betahat = H varphi.
        if (~isfield(restrict,'H'))||(size(restrict.H,2)<r)
            error('VECM_H0: Type B hypotheses needs to specify H with fewer columns than rows!');
        end
        if (size(R1t,2)~= size(restrict.H,1))
            error('VECM_H0: Type B hypotheses needs to conform with Joh_j specification of deterministics!');
        end
        H = restrict.H;
        R1H = R1t * H;
        [alphahat,phihat,~,~]= est_alpha_beta(R0t,R1H,r,Joh_j);
        betahat = H*phihat
        % restricted residuals
        res_r = R0t - R1t*betahat*alphahat';
        LL(2) =  (T-2)*log(det(res_r'*res_r/(T-2)))+(T-2)*s;
        df =  r*(s-size(H,2)+plus1);
    otherwise % must be 'C' here.
        if (~isfield(restrict,'b'))
            error('VECM_H0: Type C hypotheses needs to specify b!');
        end
        b = restrict.b;
        [Q,~]=qr(b);
        rb = size(b,2);
        b_perp = Q(:,rb+1:end);
        R1tb = R1t*b;
        R1tbp = R1t - R1tb*(R1tb\R1t);
        R0tbp = R0t - R1tb*(R1tb\R0t);

        R1tbp = R1tbp*b_perp;
        sb= size(b,2); 
        [alphahat,phihat,~,~]= est_alpha_beta(R0tbp,R1tbp,r-sb,Joh_j);
        betahat = [b,b_perp*phihat];
        % reestimate alphahat
        alphahat = ((R1t *betahat)\R0t)';

        % restricted residuals
        res_r = R0t - R1t*betahat*alphahat';
        LL(2) =  (T-2)*log(det(res_r'*res_r/(T-2)))+(T-2)*s;
        df = sb*(s-r+plus1);
end

%%%%%%%%%%%%%%%%% print out results of the test   %%%%%%%%%%%%
fprintf('\n VECM test for hypotheses on beta. \n\n');
fprintf(' T: %d, s: %d, k: %d. \n', T,s,k);

det_lab = {'trend and constant','restricted trend','constant','restricted constant','none'};

switch restrict.type
    case 'A'
        type_lab = 'H_0: beta = beta_0';
    case 'B'
        type_lab = 'H_0: beta = H phi';
    otherwise
        type_lab ='H_0: beta = [b,beta_2]';
end

fprintf('Type of restriction: %s. Type of deterministics: %s \n', type_lab, det_lab{Joh_j}); 

dLL = diff(LL);
pval = 1-chis_cdf(dLL,df);
crits =chis_inv([0.99,0.95,0.9],df);

fprintf('\n Test statistic: %3.2f, dfs: %d, p-value: %1.2f. \n',dLL,df,pval); 
fprintf('Critical values (1pct|5pct|10pct): %3.2f (1pct), %3.2f (5pct), %3.2f (10pct). \n',crits);


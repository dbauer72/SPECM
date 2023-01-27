function [LL,alphahat,betahat,df] = RH_specm_H0_Sij(R0t,R1t,r,Joh_j,rest,restrict)
% Ribarits_Hanzon idea for testing resrtictions on beta in RH_SPECM form
%
% SYNTAX: [LL,alphahat,betahat,df] = RH_specm_H0_Sij(R0t,R1t,r,Joh_j,restrict);
%
% INPUT: R0t  ... Txs concentrated left hand side.  
%        R1t  ... Tx(s+plus1) concentrated regressors.
%        r    ... rank of alpha beta'
%        Joh_j ... integer; specification of deterministics in the I(1)
%                  case of Johansen. 
%        restrict ... specifies H0.
%
% OUTPUT: LL   ... 2x1 minimized log likelihood under H0 and (H0 u H1).
%         alphahat, betahat ... sxr matrices under the restriction
%
% COMMENT: 
%          + implements estimation for three different restrictions:
%           (A) beta = restrict.beta.
%           (B) beta = restrict.H* phi
%           (C) beta = [restrict.b,beta_2]
%           restrict.type = "A"|"B"|"C".
%
% REMARK: same as RH_specm_H0, but starts after R0t, R1t are
% calculated. Thus can also be used for MFI processes. 
%
% AUTHOR: dbauer, 16.2.2022.

LL = [0,0];

plus1 = size(R1t,2)-size(R0t,2);
[T,s] = size(R0t);
% estimation without taking restrictions into account
[alphahat,betahat,~,~]= est_alpha_beta(R0t,R1t,r,Joh_j);
% unrestricted residuals
res_u = R0t - R1t*betahat*alphahat';

LL(1) =  (T)*log(det(res_u'*res_u/(T)))+(T)*s;

% estimation under H0.
switch restrict.type
    case 'A' % betahat = beta.
        if ~isfield(restrict,'beta')
            error('RH_specm_H0: Type A hypotheses needs to specify beta!');
        end
        betahat = restrict.beta;
        %r = size(betahat,2);
        R1beta = R1t * betahat;
        alphahat = (R1beta\R0t)';
        res_r = R0t - R1beta*alphahat';
        LL(2) =  (T)*log(det(res_r'*res_r/(T)))+(T)*s;
        df = r*(s+plus1-r);
    case 'B' % betahat = H varphi.
        if (~isfield(restrict,'H'))||(size(restrict.H,2)<r)
            error('RH_specm_H0: Type B hypotheses needs to specify H with fewer columns than rows!');
        end
        if (size(R1t,2)~= size(restrict.H,1))
            error('RH_specm_H0: Type B hypotheses needs to conform with Joh_j specification of deterministics!');
        end
        H = restrict.H;
        R1H = R1t * H;
        [alphahat,phihat,~,~]= est_alpha_beta(R0t,R1H,r,Joh_j);
        betahat = H*phihat;
        % restricted residuals
        res_r = R0t - R1t*betahat*alphahat';
        LL(2) =  (T)*log(det(res_r'*res_r/(T)))+(T)*s;
        df =  r*(s-size(H,2)+plus1);
    otherwise % must be 'C' here.
        if (~isfield(restrict,'b'))
            error('RH_specm_H0: Type C hypotheses needs to specify b!');
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
        LL(2) =  (T)*log(det(res_r'*res_r/(T)))+(T)*s;
        df = sb*(s-r+plus1);
end

%%%%%%%%%%%%%%%%% print out results of the test   %%%%%%%%%%%%
fprintf('\n RH_SPECM test for hypotheses on beta. \n\n');
fprintf(' T: %d, s: %d. \n', T,s);

det_lab = {'trend and constant','restricted trend','constant','restricted constant','none'};

switch restrict.type
    case 'A'
        type_lab = 'H_0: beta = beta_0.';
    case 'B'
        type_lab = 'H_0: beta = H phi.';
    otherwise
        type_lab ='H_0: beta = [b,beta_2].';
end

fprintf('Type of restriction: %s. Type of deterministics: %s \n', type_lab, det_lab{Joh_j}); 
if rest == 'y'
    fprintf('Restricted Ribarits-Hanzon SPECM is estimated.\n');
else
    fprintf('Unrestricted Ribarits-Hanzon SPECM is estimated.\n');
end

dLL = diff(LL);
pval = 1-chis_cdf(dLL,df);
crits =chis_inv([0.99,0.95,0.9],df);

fprintf('\n Test statistic: %3.2f, dfs: %d, p-value: %1.2f. \n',dLL,df,pval); 
fprintf('Critical values (1pct|5pct|10pct): %3.2f (1pct), %3.2f (5pct), %3.2f (10pct). \n',crits);


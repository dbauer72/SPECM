function [LL,alphahat,betahat,df] = RH_vecm_Ha(y,Abar,B,r,rest,Joh_j,restrict);
% Ribarits_Hanzon idea for testing resrtictions on alpha using the RH_SPECM form
%
% SYNTAX: [LL,alphahat,betahat,df] = RH_specm_Ha(y,Abar,B,r,rest,Joh_j,restrict);
%
% INPUT: y    ... Txs data matrix,
%        Abar ... nxn matrix A-BC.
%        B    ... nxs matrix
%        r    ... rank of alpha beta' 
%        rest ... 'y'|'n' restricted version.
%        Joh_j ... integer; indicating the version of Johansen's
%        deterministics.
%        restrict ... specifies H0.
%
% OUTPUT: LL   ... 2x1 minimized log likelihood under H0 and (H0 u H1).
%         alphahat, betahat ... sxr matrices under the restriction
%
% COMMENT: + The restriction C*(I-Abar)^(-1)*B = I+alpha beta' is imposed by
% rest = 'y'.
%          + implements estimation for three different restrictions:
%           (A) alpha = restrict.A psi.
%           (B) alpha = [restrict.a, restrict.a_perp phi]
%           (C) alpha = restrict.A psi, beta = restrict.H phi
%           restrict.type = "A"|"B"|"C".
%          + for Joh_j = 2,4 the additional regressor in beta is included
%          at the head not at the tail!
%
% AUTHOR: dbauer, 9.2.2022.

LL = [0,0];
% prepare regressions
if nargin< 6
    Joh_j = 5; % no deterministic terms.
end

if nargin <7
    error('RH_specm_H0: Restrictions need to be specified. See help for details.');
end

[T,s] = size(y);
n = size(Abar,1);

switch Joh_j
    case {2,4} % restricted linear trend
        plus1=1;
    otherwise
        plus1 = 0;
end

%betahat = zeros(s+plus1,r);
%alphahat = betahat;

if (n<s)&&(rest == 'y')
    message('RH_specm_H0: n<s and restricted estimation cannot both be used for testing hypothesis on beta. Changing "rest" to "n"');
    rest = 'n';
end

[R0t,R1t] = concentrate_RH(y,Abar,B,rest,Joh_j);


% estimation without taking restrictions into account
[alphahat,betahat,~,~]= est_alpha_beta(R0t,R1t,r,Joh_j);
% unrestricted residuals
res_u = R0t - R1t*betahat*alphahat';

LL(1) =  (T-2)*log(det(res_u'*res_u/(T-2)))+(T-2)*s;

% estimation under H0.
switch restrict.type
    case 'A' % hypothesis of the form H_0: alpha = A phi.
        % change hypotesis to orthonormal matrix, because handling is
        % easier.
        A = restrict.A; 
        sa = size(A,2);
        [Q,R]= qr(A);
        A = Q(:,1:sa); 
        Aperp = Q(:,sa+1:end);
        iSapap = inv(Aperp'*R0t'*R0t*Aperp);
        Rattilde = R0t*(A - Aperp*iSapap * Aperp'*R0t'*R0t*A);
        R1ttilde = R1t - R0t*Aperp*iSapap*Aperp'*R0t'*R1t;
        [xihat,betahat,~,~]= est_alpha_beta(Rattilde,R1ttilde,r,Joh_j);
        alphahat = A*xihat;

        res_r2 = [R0t*Aperp,Rattilde - R1ttilde*betahat*xihat'];
        res_r = R0t - R1t*betahat*alphahat';
        LL(2) =  (T-2)*log(det(res_r'*res_r/(T-2)))+(T-2)*s;
        LL(3) =  (T-2)*log(det(res_r2'*res_r2/(T-2)))+(T-2)*s;
        
        df = r*(s-sa);
    case 'B' % hypothesis of the form H_0: alpha = [a,tau]. 
        a = restrict.a;
        sa = size(a,2);
        [Q,R]=qr(a);
        a = Q(:,1:sa);
        aperp = Q(:,sa+1:end);
        % eq (8.19)
        Rapt = R0t*aperp;
        [psihat,beta2hat,~,~]= est_alpha_beta(Rapt,R1t,r-sa,Joh_j);
        res_r(:,(sa+1):s) = Rapt - R1t*beta2hat*psihat';
        % eq (8.20) 
        Rat = R0t*a; 
        Zt = [R1t,Rapt];
        cc = (Zt\Rat)';
        res_r(:,1:sa) = Rat - Zt*cc';
        
        % reconstruct alphahat and betahat
        alphahat = [a,aperp*psihat];
        Omegahat = res_r'*res_r/(T-2);
        omegahat = Omegahat(1:sa,(sa+1):end)*inv(Omegahat((sa+1):end,(sa+1):end));
        betahat = [cc(:,1:size(R1t,2))'+beta2hat*psihat'*omegahat',beta2hat];

        % residuals 
        res_r2 = R0t - R1t*betahat*alphahat';
        LL(2) =  (T-2)*log(det(res_r'*res_r/(T-2)))+(T-2)*s;
        LL(3) =  (T-2)*log(det(res_r2'*res_r2/(T-2)))+(T-2)*s;
        
        df = (s-r)*sa;
    otherwise % both alpha and beta are restricted
        % change hypotesis to orthonormal matrix, because handling is
        % easier.
        A = restrict.A; 
        sa = size(A,2);
        [Q,R]= qr(A);
        A = Q(:,1:sa); 
        Aperp = Q(:,sa+1:end);
        iSapap = inv(Aperp'*R0t'*R0t*Aperp);
        Rattilde = R0t*(A - Aperp*iSapap * Aperp'*R0t'*R0t*A);
        R1ttilde = R1t - R0t*Aperp*iSapap*Aperp'*R0t'*R1t;

        % this is where the restriction on beta is imposed analogous to the
        % type 'B' restrictions on beta. 
        if (~isfield(restrict,'H'))||(size(restrict.H,2)<r)
            error('RH_specm_Ha: Type C hypotheses needs to specify H with fewer columns than rows!');
        end
        if (size(R1t,2)~= size(restrict.H,1))
            error('RH_specm_Ha: Type B hypotheses needs to conform with Joh_j specification of deterministics!');
        end
        H = restrict.H;
        R1H = R1ttilde * H;
        [xihat,phihat,~,~]= est_alpha_beta(Rattilde,R1H,r,Joh_j);
        betahat = H*phihat;
        alphahat = A*xihat;

        res_r2 = [R0t*Aperp,Rattilde - R1ttilde*betahat*xihat'];
        res_r = R0t - R1t*betahat*alphahat';
        LL(2) =  (T-2)*log(det(res_r'*res_r/(T-2)))+(T-2)*s;
        LL(3) =  (T-2)*log(det(res_r2'*res_r2/(T-2)))+(T-2)*s;
        
        df =  r*(s-size(H,2)+plus1)+r*(s-sa);
end

%%%%%%%%%%%%%%%%% print out results of the test   %%%%%%%%%%%%
fprintf('\n RH_SPECM test for hypotheses on alpha. \n\n');
fprintf(' T: %d, s: %d, n: %d. \n', T,s,n);

det_lab = {'trend and constant','restricted trend','constant','restricted constant','none'};

switch restrict.type
    case 'A'
        type_lab = 'H_0: alpha = A psi.';
    case 'B'
        type_lab ='H_0: alpha = [a,alpha_2].';
    otherwise
        type_lab ='H_0: alpha = A psi, beta = H phi.';        
end

fprintf('Type of restriction: %s. Type of deterministics: %s \n', type_lab, det_lab{Joh_j}); 
if rest == 'y'
    fprintf('Restricted Ribarits-Hanzon SPECM is estimated.\n');
else
    fprintf('Unrestricted Ribarits-Hanzon SPECM is estimated.\n');
end

dLL = diff(LL(1:2));
pval = 1-chis_cdf(dLL,df);
crits =chis_inv([0.99,0.95,0.9],df);

fprintf('\n Test statistic: %3.2f, dfs: %d, p-value: %1.2f. \n',dLL,df,pval); 
fprintf('Critical values (1pct|5pct|10pct): %3.2f (1pct), %3.2f (5pct), %3.2f (10pct). \n',crits);



function [alphahat,betahat,Q,rho1,alpha_o,beta_o,Va]= est_alpha_beta(R0t,R1t,r,Joh_j);
% function estimates Pi = alphahat beta' subject to the rank restriction
% rank(Pi)=r using the Johansen approach and concentrated observations
% R_{0t}, R_{1t}. 
% 
% SYNTAX: [alphahat,betahat,Q,rho1,alpha_o,beta_o,Va]= est_alpha_beta(R0t,R1t,r,Joh_j);
%
% INPUT: R0t   ... Txs real matrix; left hand side vars of Johansen. 
%        R1t   ... Txs real matrix. right hand side vars.
%        r     ... integer; rank restriction of Pi.
%        Joh_j ... integer; indicating the deterministic terms to include
%                  according to the Johansen nomenclature: (1: unrestricted linear
%                  trend in VECM; 2: restriced linear trend; 3: unrestricted
%                  constant; 4: restricted constant; 5: no deterministics). 
%
% OUTPUT: alphahat, betahat ... s x r real matrices.
%
% COMMENTS: solves the RRR problem R0t = alpha beta' R1t + et 
% subject to rank(alpha beta')=r.
% 
% AUTHOR: dbauer, 28.11.2017.
% 

% differentiate according to Johansen determinstic specification
[T,s] = size(R0t); 

switch Joh_j
    case {2,4}
        plus1 = 1;
    otherwise
        plus1 = 0;
end;

% solve EV problem to obtain betahat.
S11 = R1t'*R1t;
S01 = R0t'*R1t;
S00 = R0t'*R0t;

% estimate alphahat
me = min(eig(S11));
if me<10^(-4)
    mm = size(S11,1);
    S11 = S11+eye(mm)*(10^(-4)-me);
end;
me = min(eig(S00));
if me<10^(-4)
    mm = size(S00,1);
    S00 = S00+eye(mm)*(10^(-4)-me);
end;

c11 = chol(S11); % c11'*c11 = S11
c00= chol(S00); % c00'*c00 = S00
B01 = inv(c00)'*S01*inv(c11);
[U,Q,V] = svd(B01);
Q = diag(Q);

% catch, if called with r==0 or r==s. 
if r==0 
    % in this case rank is restricted to zero. 
    alphahat = zeros(s,0);
    betahat = zeros(s+plus1,0);
    rho1 = zeros(0,plus1);
    %Q = 0;
    alpha_o = eye(s);
    beta_o = eye(s);
    Va = [];
    return
end
if r==s
    % in this case there are not restrictions on the rank -> unrestricted
    % estimate. 
    alphahat = eye(s); 
    betahat = (S01*inv(S11))';
    if plus1>0 
        rho1 = betahat(plus1,:);
    else
        rho1 =0;
    end;
    alpha_o = zeros(s,0);
    beta_o = zeros(s,0);
    Va = [];
    return
end


% plus1=0;
% switch Joh_j
%     case 1
%         %yh=[1:T]'/T;
%         %yh=yh.^2;
%         %R1t = [yh,R1t];
%         %plus1=1;
%         R1t = detrend(R1t,'linear');
%         R0t = detrend(R0t,'linear');
%     case 2 % restricted linear trend
%         yh=[1:T]'/T;
%         R1t = [yh,R1t];
%         plus1=1;
%         R1t = detrend(R1t,0);
%         R0t = detrend(R0t,0);
%     case 3 % unrestricted constant
%         %yh=ones(T-1,1);
%         %R1t = [yh,R1t];
%         R1t = detrend(R1t,0);
%         R0t = detrend(R0t,0);
%     case 4 % restricted constant.
%         yh=ones(T,1);
%         R1t = [yh,R1t];
%         plus1=1;
%     otherwise        
% end


alphahat = c00'*U(:,1:r)*diag(Q(1:r));
betahat = inv(c11)*V(:,1:r);

beta_o = S11*inv(c11)*V(:,(r+1):end);
alpha_o = inv(S00)*S01*inv(c11)*V(:,(r+1):end);

% normalize beta to orthonormal block column. 
[~,R]= qr(betahat);
Tr = R(1:r,1:r)';
alphahat = alphahat*Tr;
betahat = betahat*inv(Tr)';

switch Joh_j %additional vars included at first position due to deterministics? 
    case {2,4}
        rho1 = betahat(1,:);
        %betahat = betahat(2:end,:);
    otherwise
        rho1 = zeros(r,0);
end

% estimate variance of alpha for inference
Vbb = betahat'*S11*betahat;
Omega= S00 - alphahat*Vbb*alphahat';
Va = kron(Omega,inv(Vbb));

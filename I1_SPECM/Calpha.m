function [alphahat,Valpha] = Calpha(y,Abar,B,betahat,Joh_j,correct);
% Re-estimation of alpha after estimating beta using C(alpha) method.
%
% SYNTAX: [alphahat] = Calpha(y,Abar,B,betahat,Joh_j);
%
% INPUT: y    ... Txs data matrix,
%        Abar ... nxn matrix A-BC.
%        B    ... nxs matrix
%        betahat ... (s+plus1)xr matrix of estimated cointegrating relations.
%        Joh_j ... integer; indicating the version of Johansen's
%        deterministics. 
%        correct ... indicator, if Calpha correction should be applied.
%
% OUTPUT: alphahat ... sxr matrices.
%
% COMMENT: C(alpha) means that the estimates (Abar,B) are corrected for. 
%
% AUTHOR: dbauer, 10.12.2021. 

% prepare regressions 
if nargin< 6
    Joh_j = 5; % no deterministic terms. 
end;

[T,s] = size(y);
n = size(Abar,1);
% renormalize for calculation of derivatives
Controll = B;
if n>s
    while (size(Controll,2)<n)
        Controll = [B,Abar*Controll];
    end
end

Trafo = Controll(:,1:n); 
Abar = inv(Trafo)*Abar*Trafo;
B = inv(Trafo)*B; 

% steps: calculate tilde x_t: filtered version of Delta y.
Z0t = diff(y);
%Z0t = detrend(Z0t,0);
Z1t = y(1:end-1,:);

Z2t = zeros(T,n);

tB = B';
tA = Abar';
for t=2:T
    Z2t(t,:)=Z2t(t-1,:)*tA + Z0t(t-1,:)*tB;
end;
%Z2t = ltitr(Abar,B,Z0t,zeros(n,1));

Z2t = Z2t(2:end-1,:);
Z0t = Z0t(2:end,:);
Z1t = Z1t(2:end,:);

plus1 = 0;
switch Joh_j
    case 1
        Z0t = detrend(Z0t,1);
        Z1t = detrend(Z1t,1);
        Z2t = detrend(Z2t,1);
    case 2 % restricted linear trend
        Z0t = detrend(Z0t,0);
        Z1t = [[1:T-2]',Z1t];
        Z1t = detrend(Z1t,0);
        Z2t = detrend(Z2t,0);
        plus1=1;
    case 3 % unrestricted constant
        Z0t = detrend(Z0t,0);
        Z1t = detrend(Z1t,0);
        Z2t = detrend(Z2t,0);
    case 4
        Z1t = [ones(T-2,1),Z1t];
        plus1=1;
end

% check dimensions
if size(betahat,1)~= size(Z1t,2)
    if size(betahat,1)== size(Z1t,2)-plus1
        betahat = [betahat;zeros(1,size(betahat,2))];
        warning("Betahat does not include additional deterministic component, adding zeros.")
    else
        error('Wrong dimension of betahat!');
    end
end

% set up regressors
Xt = [Z1t*betahat,Z2t];
r = size(betahat,2); 

% first regression in order to obtain Chat. 
alphaC = (Xt\ Z0t)'; 
Chat = alphaC(:,(r+1):end);

% generate additional derivatives 
dZ2t = gen_dZ2t(Z0t,Z2t,Abar,B);

% iterate over output components 
res = 0*Z0t;
Xts_f = zeros(size(res,1),0);
for j = 1:s
    Xts= Xt;
    if correct==1
        for a=1:n
            Xts = [Xts,dZ2t(:,(a-1)*n+[1:n])*Chat(j,:)'];
        end
    end
    alphaC_dC(j,:) = (Xts\ Z0t(:,j))'; 
    res(:,j) = Z0t(:,j) - Xts*alphaC_dC(j,:)';
    Xts_f = [Xts_f,Xts*inv(Xts'*Xts)];
end

% dimensions
nts = size(alphaC_dC,2);
Valpha = zeros(r*s,r*s);
Omega = res'*res/size(res,1);

for a=1:s
    for b=1:s
        Valpha(r*(a-1)+[1:r],r*(b-1)+[1:r]) = Omega(a,b)*Xts_f(:,r*(a-1)+[1:r])'*Xts_f(:,r*(b-1)+[1:r]);
    end
end

alphahat = alphaC_dC(:,1:r); 



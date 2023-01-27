function [LL,alphahat,betahat,res,S,th,Va] = VECM_I1(y,k,r,Joh_j);
% VECM form estimation of AR systems in the I(1) case following Johansen's
% VECM formulation.
%
% SYNTAX: [LL,alphahat,betahat,res] = VECM(y,k,r,Joh_j);
%
% INPUT: y    ... Txs data matrix,
%        k    ... lag order of VECM Delta terms.
%        r    ... int; coint-rank.
%        Joh_j ... deterministics of Johansen. 
%
% OUTPUT: LL   ... minimized log likelihood.
%         alphahat, betahat ... sxr matrices.
%         res ... Txs residuals
%
% REMARK: Joh_j is an integer containing the specification of the
% characteristics according to Johansen (1997). 
%  Joh_j:    1 ... no restriction. 
% AUTHOR: dbauer, 19.12.2019. 

% prepare regressions 
if nargin< 3
    Joh_j = 5; % no deterministic terms. 
end;
    
[T,s] = size(y);

Z0t = diff(y);
Z1t = y(1:end-1,:);

Z2t = zeros(T-1,s*k);

if k<0 % negative lag length implies AIC estimation of order.
    k = aicest(y,s,abs(k));
end


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

[alphahat,betahat,S,rho1,~,~,Va]= est_alpha_beta(R0t,R1t,r,Joh_j);

% residuals
Z0t = diff(y);
%Z1t = y(1:end-1,:);
res = Z0t(k+1:end,:) - Z1to*betahat*alphahat';

switch Joh_j
    case 1 % no restriction, full constant and linear trend 
        
    %    res = detrend(res,'linear');
    case 2 % restricted linear trend
        %yh=[1:(T-k-1)]'/(T-k-1);
       % res = res - Z1to(:,plus1)*rho1*alphahat';
    %    res = detrend(res,0);
    %case 3 % unrestricted constant
    %    res = detrend(res,0);
    case 4 % restricted constant.
        %yh=ones(T-k-1,1);
     %   res = res - Z1to(:,plus1)*rho1*alphahat';
    %otherwise
end
 
Zst = [Z2to,dt(k+2:end,:)];
Phi = (Zst\res)';
res = res - Zst*Phi';

% convert VECM to VAR formulation. 
A = zeros(s,(k+2)*s);
A(:,1:s)=eye(s);
A(:,s+1:2*s) = -alphahat*betahat(1:s,:)'-eye(s)-Phi(:,[1:s]);
for j=2:k
    A(:,j*s+[1:s])= -Phi(:,(j-1)*s+[1:s])+Phi(:,(j-2)*s+[1:s]);
end

A(:,(k+1)*s+[1:s])=  Phi(:,(k-1)*s+[1:s]);

% fill in values of th structure
th = theta_urs();
th.a = A;
th.b = eye(s);
% differentiate deterministic cases from Johansen. 

switch Joh_j
    case 1 % no restrictions 
        B = Phi(:,end-1:end);
        m= 2;
    case 2 % restricted linear trend        
        %res = res - Z1t(:,plus1)*rho1*alphahat';
        B(:,1)= Phi(:,end);
        if size(alphahat,2)>0
            B(:,2)= alphahat*rho1';
        else 
            B(:,2)=0;
        end;
        m=2;
    case 3 % unrestricted constant
    	B(:,1)= Phi(:,end);
        m=1;
    case 4 % restricted constant.        
        %res = res - Z1t(:,plus1)*rho1*alphahat';
        if size(alphahat,2)>0
            B= alphahat*rho1';
        else 
            B=zeros(s,1);
        end;
        m=1;
    otherwise
        B = zeros(s,0);
        m=0;
end


th.Omega = res'*res/(T-k-1);
th.which = 'poly'; 
th.m = m;
th.d = B;

% mark as unit root process
th.ur = 'I(1)';
th.urs = s-r; 

% log likelihood
LL =  (T-k-1)*log(det(res'*res/(T-k-1)))+(T-k-1)*s;

res = [NaN(k+1,s);res];


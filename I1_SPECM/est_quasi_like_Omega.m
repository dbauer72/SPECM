function [Omegah,qlike,tres] = est_quasi_like_Omega(th,y,s,m);
% estimates the innovation variance from prediction error estimates. 
%  quasi likelihood for theta structure theta.
%
% SYNTAX: [Omegah,qlike,tres] = est_quasi_like_Omega(th,y,s,m);
%
% INPUT:    th ... theta structure of estimated system
%           y  ... Tx(s+m) matrix of observations.
%           s  ... integer; dimension of endogenous vars.
%           m  ... integer; dimension of exogenous vars. 
%
% OUTPUT:   Omegah ... sxs matrix of estimated innovation variance
%           qlike ... real; -2/T log of Gaussian likelihood. 
%           tres  ... Txs matrix of residuals.
%
% REMARK: Omegah = tres'*tres/T.
%
% AUTHOR: dbauer, 19.12.2019. 



qlike = 0;

Omega = th.Omega;

tilA = th.A;
tilK = th.K;
tilC = th.C;
D = th.D; 

n = size(tilA,1);
%
% convert params to matrices
dt = y(:,s+1:end);
m = size(dt,2);
y = y(:,1:s);

T = size(y,1);
ty = y - dt*D';


% calculate P_bull
Q = tilK*Omega*tilK';
Q = (Q+Q')/2;
P0 = 0*eye(n);

% --- initialize the Kalman filter ---
x0= zeros(n,1);
xf = x0 + P0*tilC'*inv(tilC*P0*tilC'+Omega)*(y(1,:)'-tilC*x0);
P0g0 = P0-P0*tilC'*inv(tilC*P0*tilC'+Omega)*tilC*P0;

% --- run the Kalman filter ---

xf = tilA*xf;
Pkg1 = tilA*P0g0*tilA' + Q;
tres = ty*0;
tres(1,:)=ty(1,:);
Omegat= tilC*P0*tilC'+Omega;
qlike = qlike +  log(det(Omegat)) + tres(1,:)*inv(Omegat)*tres(1,:)'; 

for t=2:T % filter equations 
    Omegat = tilC*Pkg1*tilC'+Omega;
    iOm = inv(Omegat);
    tres(t,:)= ty(t,:) - xf'*tilC';
    Kt = (tilA*Pkg1*tilC'+tilK*Omega)*iOm;
    xf = tilA*xf+ Kt*tres(t,:)';
    Pkg1 = (tilA * Pkg1 *tilA') +Q - (Kt*Omegat*Kt');    
    Pkg1 = (Pkg1+Pkg1')/2;
    % update likelihood 
    qlike = qlike + log(det(Omegat)) + tres(t,:)*inv(Omegat)*tres(t,:)'; 
end

%qlike = qlike/T;
%Omegah = tres'*tres/T;

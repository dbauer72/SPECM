function [qlike,tres] = cal_quasi_like_th(param,th,y,s,m,n,c,Pbullet,restrict);
% calculates the quasi likelihood for parameter vectors param using Omega
% and D from th.
%
% SYNTAX: [qlike,tres] = cal_quasi_like_th(param,th,y,s,m,n,c,Pbullet,restrict);
%
% INPUTS:  param ... d x 1 vector of parameter values (fed into param2th,
%                       see there for description)
%          th    ... theta structure for state space system;
%          y     ... T x s matrix of observations
%          s     ... integer; dimension of endogenous vars
%          m     ... integer; dimension of exogenous vars (e.g.
%                     deterministics)
%          n     ... integer; state dimension
%          c     ... integer; number of common trends.
%          Pbullet ... indicator var; if >0: state started with stationary
%                   variance; if 0: state initialized with zero (corresponds to prediction error).
%          restrict ... structure containing restrictions to be passed on
%                       to param2syst.
%
% OUTPUT:   qlike ... real; -2/T log Gaussian likelihood.
%           tres  ... Txs; matrix of residuals.
%
% REMARKS:
%         + restrict.Omega: if this field is included, Omega is fixed not
%                         included in parameter vector.
%         + c>0: In this case the observations are transformed to stationary
%         values for t>1.
%         + large penalization of eigenvalues larger than 0.99 of tilde A (unit
%         roots) and tilde A-KC (minimum phase).
%
% AUTHOR: dbauer, 27.10.2019.

if nargin<8
    restrict.det_res =0;
    restrict.scale = ones(length(param),1);
end;

if nargin<7 % Pbullet: start with x_1 = x_bullet?
    Pbullet = 0;
end;

qlike = 0;
% extract D and Omega from th
D = th.D;
Omega = th.Omega;

% calculate system from parameters.
the = param2th(param,s,n,c,restrict);

A = the.A;
K = the.K;
C = the.C;

dt = y(:,s+1:end);
m = size(dt,2);
y = y(:,1:s);

T = size(y,1);
y = y - dt*D';

if c>0 % there is an integration going on!
    % find tilde k. via C1:
    C1 = C(:,1:c);
    Pi = C1*inv(C1'*C1)*C1';
    Cbull = [C(:,c+1:end)];
    Kbull= K(c+1:end,:);
    K1 = K(1:c,:);
    
    % calculate tilde y_t = y_t - Pi y_{t-1}
    ty = y;
    ty(2:end,:)=ty(2:end,:)-y(1:end-1,:)*Pi;
    
    % transformed matrices
    tilA = [zeros(c,c),-C1'*Cbull;A(c+1:end,:)];
    tilK = [K1-C1';Kbull];
    tilC = C;
else % no integration
    tilA = A;
    tilK= K;
    tilC=C;
    ty = y;
end;

if Pbullet>0 % if P_bull to be included: system cannot be unstable!
    maA = max(abs(eig(tilA)));
    if maA >0.99
        tilA = tilA*0.99/maA;
        qlike = qlike + 10^6*(exp(maA)-exp(0.99));
    end
end
% correct for not-invertible systems: add large penalty.
makcA = max(abs(eig(tilA-tilK*tilC)));
if makcA>0.99
    tilA = tilA*0.99/makcA;
    tilK = tilK*0.99/makcA;
    qlike = qlike + 10^6*(exp(makcA)-exp(0.99));
end


% calculate P_bull
Q = tilK*Omega*tilK';
Q = (Q+Q')/2;
Qbull = Q(c+1:end,c+1:end);
if Pbullet>0
    pb = inv( eye((n-c)^2) - kron(tilA(c+1:end,c+1:end),tilA(c+1:end,c+1:end)))*Qbull(:);
    Pbull = reshape(pb,n-c,n-c);
    P0 = [0*eye(c),zeros(c,n-c);zeros(n-c,c),Pbull];
    
    P0 = (P0+P0')/2;
    % --- initialize the Kalman filter ---
    x0= zeros(n,1);
    tres = ty*0;
    tres(1,:)=ty(1,:);
    Omegat= tilC*P0*tilC'+Omega;
    Kt = (tilA*P0*tilC'+tilK*Omegat)*inv(Omegat);
    xf = tilA*x0 + Kt*tres(1,:)';
    Pkg1 = tilA*P0*tilA' + Q- Kt*Omegat*Kt'; % P(2|1)
    
    %xf = x0 + P0*tilC'*inv(tilC*P0*tilC'+Omega)*(y(1,:)'-tilC*x0);
    %P0g0 = P0-P0*tilC'*inv(tilC*P0*tilC'+Omega)*tilC*P0;
    %
    %% --- run the Kalman filter ---
    %
    %xf = tilA*xf;
    %Pkg1 = tilA*P0g0*tilA' + Q;
    %tres = ty*0;
    %tres(1,:)=ty(1,:);
    %Omegat= tilC*P0*tilC'+Omega;
    
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
else
    P0 = 0*eye(n);
    % --- initialize the Kalman filter ---
    x0= zeros(n,1); % x(1|0)
    Abar = tilA-tilK*tilC;
    % --- run the filter for the inverse  transfer function ---
    Ts = size(ty,1);
    xt = ltitr(Abar,tilK,ty(1:Ts,:),x0);
    tres = ty(1:(Ts),:)-xt*tilC';
    [tres,~] = est_initial_val(tres,Abar,tilK,tilC);
    Omegat =  tres(1:(Ts),:)'*tres(1:(Ts),:)/(Ts);
    
    % update likelihood
    qlike = qlike + Ts*log(det(Omegat)) + s*Ts;
end
%qlike = qlike/T;

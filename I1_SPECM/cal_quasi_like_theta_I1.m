function [qlike,tres] = cal_quasi_like_theta_I1(th,y,s,m,n,c,Pbullet,restrict);
% calculates the quasi likelihood for theta structure theta.
%
% SYNTAX: [qlike,tres] = cal_quasi_like_theta_I1(th,y,s,m,n,c,Pbullet,restrict);
%
% INPUTS:  th ... theta structure of system to evaluate the quasi
%                    likelihood.
%          y     ... T x s matrix of observations
%          s     ... integer; dimension of endogenous vars
%          m     ... integer; dimension of exogenous vars (e.g.
%                     deterministics)
%          n     ... integer; state dimension
%          c     ... integer; number of common trends.
%          Pbullet ... indicator var; if >0: state started with stationary
%                   variance; if 0: state initialized with zero (corresponds to prediction error).
%          restrict ... not used; included for consistency in inputs for
%                    function.
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


if nargin<7 % Pbullet: start with x_1 = x_bullet?
    Pbullet = 0;
end;

qlike = 0;
% extract matrices from theta
A=th.A;
K=th.K;
C = th.C;
D = th.D;
Omega = th.Omega;

% extract parameters for Omega
%parom = param(1:(s*(s+1)/2));
%param(1:(s*(s+1)/2))=[];
%
%Omega = fill_lowtri(parom,s);
%
% convert params to matrices
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
    if maA >0.999
        tilA = tilA*0.999/maA;
        qlike = qlike + 10^6*(exp(maA)-exp(0.999));
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
    
    
    P0 = (P0+P0')/2; % P(1|0).
    % --- initialize the Kalman filter ---
    x0= zeros(n,1); % x(1|0)
    
    % --- run the Kalman filter ---
    tres = ty*0;
    tres(1,:)=ty(1,:); % e(1)
    Omegat= tilC*P0*tilC'+Omega; % Omega(1|0)
    Kt = (tilA*P0*tilC'+tilK*Omegat)*inv(Omegat);
    xf = tilA*x0 + Kt*tres(1,:)';
    Pkg1 = tilA*P0*tilA' + Q- Kt*Omegat*Kt'; % P(2|1)
    
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
    tres = ty(1:Ts,:)-xt*tilC';

    [tres,~] = est_initial_val(tres,Abar,tilK,tilC);
    Omegat =  tres(1:Ts,:)'*tres(1:Ts,:)/(Ts);
    
    % update likelihood
    qlike = qlike + Ts*log(det(Omegat)) + s*Ts;
end
%qlike = qlike/T;



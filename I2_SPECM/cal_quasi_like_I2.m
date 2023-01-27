function [qlike,tres] = cal_quasi_like_I2(param,y,s,m,n,ds,Pbullet,restrict);
% calculates the quasi likelihood for theta structure theta.
%
% SYNTAX: [qlike,tres] = cal_quasi_like(param,y,s,m,n,c,Pbullet,restrict);
%
% INPUTS:  param ... d x 1 vector of parameter values (fed into param2syst,
%                       see there for description)
%          y     ... T x s matrix of observations
%          s     ... integer; dimension of endogenous vars
%          m     ... integer; dimension of exogenous vars (e.g.
%                     deterministics)
%          n     ... integer; state dimension
%          ds     ... pair of integers; number of [I(2),I(1)] common trends.
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

if isfield(restrict,'Omega')
    Om = restrict.Omega;
    paraomi = extr_lowtri(Om);
    param = [paraomi(:);param];
end;

qlike = 0;
% extract matrices from parameters 
[A,K,C,D,Omega] = param2syst_I2(param,s,n,m,ds,restrict); 

% convert params to matrices
dt = y(:,s+1:end);
m = size(dt,2);
y = y(:,1:s);

T = size(y,1);
y = y - dt*D';

% count number of common trends
c = sum(ds)+ds(1);

if c>0 % there is an integration going on! 
    % 1st STEP: take out I(2) component.    
    % find tilde k. via C1: 
    Cbull = C(:,c+1:end);
    C1 = C(:,1:ds(1));
    C2 = C(:,ds(1)+[1:ds(1)]);
    C3 = C(:,2*ds(1)+[1:ds(2)]);

    Pi = C1*inv(C1'*C1)*C1';   
    K1 = K(1:ds(1),:);
    K2 = K(ds(1)+[1:ds(1)],:);
    K3 = K(2*ds(1)+[1:ds(2)],:);
    Kbull = K(c+1:end,:);
    
    % calculate tilde y_t = y_t - Pi y_{t-1}
    ty = y; 
    ty(2:end,:)=ty(2:end,:)-y(1:end-1,:)*Pi;

    % transformed matrices 
    sds = sum(ds);
    tilA = [eye(sds),zeros(sds,n-sds);zeros(n-sds,sds+ds(1)),[-C1'*Cbull;A(c+1:end,c+1:end)]];
    tilK = [K2;K3;-K1+K2+C1';Kbull];
    tilC = [C1+C2,C3,-C1,Cbull];
    
    % 2nd STEP: take out I(1) components from here.
    Ctil = [C1+C2,C3];
    Pi123 = Ctil*inv(Ctil'*Ctil)*Ctil';
    
    tty = ty;
    tty(2:end,:)=tty(2:end,:) - tty(1:end-1,:)*Pi123;
    
    % transform matrices 
    Cbull = tilC(:,sds+1:end);
    Kbull = tilK(sds+1:end,:);
    C1n = tilC(:,1:sds);
    K1n = tilK(1:sds,:);
    
    ttilA = [zeros(sds,sds),-Ctil'*Cbull;tilA(sds+1:end,:)];
    ttilK = [K1n-C1n';Kbull];
    ttilC = tilC;
    
    
else % no integration
    ttilA = A;
    ttilK= K;
    ttilC=C;
    tty = y;
end;

if Pbullet>0 % if P_bull to be included: system cannot be unstable!
    maA = max(abs(eig(ttilA)));
    if maA >0.99
        ttilA = ttilA*0.99/maA;
        qlike = qlike + 10^6*(exp(maA)-exp(0.99));
    end
end
% correct for not-invertible systems: add large penalty.
makcA = max(abs(eig(ttilA-ttilK*ttilC)));
if makcA>0.99
    ttilA = ttilA*0.99/makcA;
    ttilK = ttilK*0.99/makcA;
    qlike = qlike + 10^6*(exp(makcA)-exp(0.99));
end

% calculate P_bull
Q = ttilK*Omega*ttilK';
Q = (Q+Q')/2;
Qbull = Q(c+1:end,c+1:end);
if Pbullet>0  
    pb = inv( eye((n-c)^2) - kron(ttilA(c+1:end,c+1:end),ttilA(c+1:end,c+1:end)))*Qbull(:);
    Pbull = reshape(pb,n-c,n-c);
    P0 = [0*eye(c),zeros(c,n-c);zeros(n-c,c),Pbull];
else 
    P0 = 0*eye(n);
end

P0 = (P0+P0')/2;
% --- initialize the Kalman filter ---
x0= zeros(n,1);
tres = tty*0;
tres(1,:)=tty(1,:);
Omegat= ttilC*P0*ttilC'+Omega;
Kt = (ttilA*P0*ttilC'+ttilK*Omegat)*inv(Omegat);
xf = ttilA*x0 + Kt*tres(1,:)';
Pkg1 = ttilA*P0*ttilA' + Q- Kt*Omegat*Kt'; % P(2|1)

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
    Omegat = ttilC*Pkg1*ttilC'+Omega;
    iOm = inv(Omegat);
    tres(t,:)= tty(t,:) - xf'*ttilC';
    Kt = (ttilA*Pkg1*ttilC'+ttilK*Omega)*iOm;
    xf = ttilA*xf+ Kt*tres(t,:)';
    Pkg1 = (ttilA * Pkg1 *ttilA') +Q - (Kt*Omegat*Kt');    
    Pkg1 = (Pkg1+Pkg1')/2;
    % update likelihood 
    qlike = qlike + log(det(Omegat)) + tres(t,:)*inv(Omegat)*tres(t,:)'; 
end

%qlike = qlike/T;

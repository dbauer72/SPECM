function [qlike,grad_q,tres,psi] = cal_quasi_like_grad(param,y,s,m,n,c,Pbullet,restrict);
% calculates the quasi likelihood for the parameter vector param including
% its gradient. 
%
% SYNTAX: [qlike,grad_q,tres,psi] = cal_quasi_like(param,y,s,m,n,c,Pbullet,restrict);
%
% INPUTS:  param ... d x 1 vector of parameter values (fed into param2syst,
%                       see there for description)
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
%           grad ...  dx1; vector of gradients. 
%           tres  ... Txs; matrix of residuals. 
%           psi   ... Txsxd; array of derivatives of innovation estimates
%           w.r.t. parameters. 
%
% REMARKS: 
%         + restrict.Omega: if this field is included, Omega is fixed not
%                         included in parameter vector.
%         + c>0: In this case the observations are transformed to stationary
%         values for t>1. 
%         + large penalization of eigenvalues larger than 0.99 of tilde A (unit
%         roots) and tilde A-KC (minimum phase). 
%         + numerically slow, only to be used for final estimated system.
%         
% AUTHOR: dbauer, 14.1.2020. 

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
np = length(param);
grad_q = zeros(np,1);
% extract matrices from parameters 

[A,K,C,D,Omega,~,gr] = param2syst(param,s,n,m,c,restrict); 

% convert params to matrices
dt = y(:,s+1:end);
m = size(dt,2);
y = y(:,1:s);

T = size(y,1);
y = y - dt*D';
for j=1:np
    dy(:,:,j)= -dt*gr(j).D';
end;

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

    for j=1:np
        barC = C1*inv(C1'*C1);
        IPi= eye(s)-Pi;        
        dPi = IPi*gr(j).C(:,1:c)*barC';
        dPi = dPi+ dPi';
        dy(2:end,:,j)= dy(2:end,:,j) - dy(1:end-1,:,j)*Pi - y(1:end-1,:)*dPi;
    end
    % transformed matrices 
    tilA = [zeros(c,c),-C1'*Cbull;A(c+1:end,:)];
    tilK = [K1-C1';Kbull];
    tilC = C;
    
    % transform derivatives 
    for j=1:np
        gr(j).A = [zeros(c,c),-gr(j).C(:,1:c)'*Cbull-C1'*gr(j).C(:,c+1:end);gr(j).A(c+1:end,:)];
        gr(j).K = gr(j).K - [gr(j).C(:,1:c)';zeros(n-c,s)];        
    end        
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
for j=1:np
    dQ(:,:,j) = gr(j).K*Omega*tilK' + tilK*gr(j).Omega*tilK' + tilK*Omega*gr(j).K';
end

if Pbullet>0  
    pb = inv( eye((n-c)^2) - kron(tilA(c+1:end,c+1:end),tilA(c+1:end,c+1:end)))*Qbull(:);
    Pbull = reshape(pb,n-c,n-c);
    P0 = [0*eye(c),zeros(c,n-c);zeros(n-c,c),Pbull];
else 
    P0 = 0*eye(n);
end

P0 = (P0+P0')/2;
% --- initialize the Kalman filter ---
x0= zeros(n,1);
xf(:,1) = x0 + P0*tilC'*inv(tilC*P0*tilC'+Omega)*(y(1,:)'-tilC*x0);
P0g0 = P0-P0*tilC'*inv(tilC*P0*tilC'+Omega)*tilC*P0;

% gradient
for j=1:np
    iP = inv(tilC*P0*tilC'+Omega);
    dOmegat(:,:,j) = gr(j).C*P0*tilC' + tilC*P0*gr(j).C' + gr(j).Omega;
    dCiP= (P0*gr(j).C'*iP-P0*tilC'*iP*dOmegat(:,:,j)*iP);
    dxf(:,1,j)= dCiP*(y(1,:)'-tilC*x0)-P0*tilC'*iP*gr(j).C*x0;
    dPg(:,:,j)= -dCiP*tilC*P0 - P0*tilC'*iP*gr(j).C*P0;
end

% --- initialize the Kalman filter ---
x0= zeros(n,1);
tres = ty*0;
tres(1,:)=ty(1,:);
Omegat= tilC*P0*tilC'+Omega;
Kt = (tilA*P0*tilC'+tilK*Omegat)*inv(Omegat);
xf(:,2) = tilA*xf(:,1) + Kt*tres(1,:)';
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


% gradient
for j=1:np
    iOm = inv(Omegat);
    dA= gr(j).A;    
    dxf(:,2,j) = tilA*dxf(:,1,j) + gr(j).A*xf(:,1);
    dPg1(:,:,j)= dA*P0g0*tilA' + tilA*dPg(:,:,j)*tilA' + tilA*P0g0*dA' + dQ(:,:,j);
    psi(1,:,j)=dy(1,:,j);
    grad_q(j)=grad_q(j)+ trace(iOm*dOmegat(:,:,j))-tres(1,:)*iOm*dOmegat(:,:,j)*iOm*tres(1,:)'+  psi(1,:,j)*iOm*tres(1,:)'  + tres(1,:)*iOm*psi(1,:,j)' ;
end

for t=2:T % filter equations 
    Omegat = tilC*Pkg1*tilC'+Omega;
    iOm = inv(Omegat);
    tres(t,:)= ty(t,:) - xf(:,t)'*tilC';
    Kt = (tilA*Pkg1*tilC'+tilK*Omega)*iOm;
    xf(:,t+1) = tilA*xf(:,t)+ Kt*tres(t,:)';
    Pkg1_old = Pkg1;
    Pkg1 = (tilA * Pkg1 *tilA') +Q - (Kt*Omegat*Kt');    
    Pkg1 = (Pkg1+Pkg1')/2;
    % update likelihood 
    qlike = qlike + log(det(Omegat)) + tres(t,:)*inv(Omegat)*tres(t,:)'; 
    
    % gradient 
    for j=1:np
        dA= gr(j).A;dC= gr(j).C; 
        dOmegat(:,:,j) = dC*Pkg1_old*tilC' + tilC*dPg1(:,:,j)*tilC'+tilC*Pkg1_old*dC' + gr(j).Omega;
        diOm = -iOm*dOmegat(:,:,j)*iOm;
        psi(t,:,j)= dy(t,:,j)- dxf(:,t,j)'*tilC' - xf(:,t)'*dC';
        dKt = (dA*Pkg1_old*tilC'+tilA*dPg1(:,:,j)*tilC'+tilA*Pkg1_old*dC'+gr(j).K*Omega+tilK*gr(j).Omega)*iOm + (tilA*Pkg1_old*tilC'+tilK*Omega)*diOm;
        dxf(:,t+1,j)= dA*xf(:,t)+ tilA*dxf(:,t,j)+ dKt*tres(t,:)'+ Kt*psi(t,:,j)';
        dPg1(:,:,j)= dA*Pkg1_old * tilA' + tilA*dPg1(:,:,j)*tilA' + tilA*Pkg1_old*dA' +dQ(:,:,j) - (dKt*Omegat*Kt' + Kt*dOmegat(:,:,j)*Kt' + Kt*Omegat*dKt');
        grad_q(j) = grad_q(j) + trace(iOm*dOmegat(:,:,j)) + psi(t,:,j)*iOm*tres(t,:)' + tres(t,:)*diOm*tres(t,:)' + tres(t,:)*iOm*psi(t,:,j)';
    end
end

%qlike = qlike/T;
%grad_q = grad_q/T;

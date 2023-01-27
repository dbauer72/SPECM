function [qlike,res] = cal_quasi_like_MFI1(param,y,s,S,n,urs,Pbullet,restrict);
% calculates the quasi likelihood for theta structure theta.
%
% SYNTAX: [qlike,res] = cal_quasi_like_theta_MFI1(th,y,s,m,n,c,Pbullet,restrict);
%
% INPUTS:  param ... parameter vector
%          y     ... T x s matrix of observations
%          s     ... integer; dimension of endogenous vars
%          S     ... integer; number of seasons.
%          n     ... integer; state dimension
%          urs   ... unit root structure.
%          Pbullet ... indicator var; if >0: state started with stationary
%                   variance; if 0: state initialized with zero (corresponds to prediction error). 
%          restrict ... not used; included for consistency in inputs for
%                    function. 
%
% OUTPUT:   qlike ... real; -2/T log Gaussian likelihood.
%           res  ... Txs; matrix of residuals. 
%
% REMARKS: 
%         + restrict.Omega: if this field is included, Omega is fixed not
%                         included in parameter vector.
%         + observations are transformed using the I(1) approach in the VS representation.  
%         + large penalization of eigenvalues larger than 0.99 of tilde A (unit
%         roots) and tilde A-KC (minimum phase). 
%         
% AUTHOR: dbauer, 14.10.2020. 


if nargin<7 % Pbullet: start with x_1 = x_bullet? 
    Pbullet = 0;
end;

qlike = 0;
% extract matrices from parameter
[A,K,C,D,Omega] = param2syst_MFI1(param,s,n,urs,restrict); 


% convert params to matrices
dt = y(:,s+1:end);
m = size(dt,2);
y = y(:,1:s);

% convert y to VS representation. 
T = size(y,1);
y = y - dt*D';

Tau = floor(T/S);
Y = zeros(Tau,s*S);
As = A^S;
Cs = zeros(s*S,n);
Ks = zeros(n,s*S);
Es = zeros(s*S,s*S);

for j=1:S
   Y(1:Tau,(j-1)*s+[1:s])=y(j:S:S*Tau,:);
   Cs((j-1)*s+[1:s],:) = C*(A^(j-1));
   Ks(:,(j-1)*s+[1:s])=(A^(S-j))*K;
   for k=1:(j-1)
       Es((j-1)*s+[1:s],(k-1)*s+[1:s])= C*(A^(j-k-1))*K;
   end
   Es((j-1)*s+[1:s],(j-1)*s+[1:s])=eye(s);
end;

Ks = Ks*inv(Es);
% add last year, if not multiple. 
if (T>Tau*S)
    Y(Tau+1,:)=NaN;
    for t=(Tau*S+1):T
        j = t-Tau*S;
        Y(Tau+1,(j-1)*s+[1:s])=y(t,:);
    end;
end

Omegas = Es*kron(eye(S),Omega)*Es';

c = sum((urs(:,2)+1).*urs(:,3));




if c>0 % there is an integration going on! 
    % find tilde k. via C1: 
    C1 = Cs(:,1:c);
    [Q,R]= qr(C1);
    C1 = Q(:,1:c);
    Pi = C1*C1';
    Cbull = [Cs(:,c+1:end)];
    Kbull= Ks(c+1:end,:);
    K1 = R(1:c,1:c)*Ks(1:c,:);

    % calculate tilde y_t = y_t - Pi y_{t-1}
    ty = Y; 
    ty(2:end,:)=ty(2:end,:)-Y(1:end-1,:)*Pi;

    % transformed matrices 
    tilA = [zeros(c,c),-C1'*Cbull;As(c+1:end,:)];
    tilK = [K1-C1';Kbull];
    tilC = [C1,Cbull];
else % no integration
    tilA = As;
    tilK= Ks;
    tilC=Cs;
    ty = Y;
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
    if (makcA>1.5)
        Omegas = eye(size(Omegas,1));
    end
    qlike = qlike + min(10^9,10^6*(exp(makcA)-exp(0.99)));
end


% calculate P_bull
Q = tilK*Omegas*tilK';
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
    Omegat= tilC*P0*tilC'+Omegas; % Omega(1|0)
    Kt = (tilA*P0*tilC'+tilK*Omegat)*inv(Omegat);
    xf = tilA*x0 + Kt*tres(1,:)';
    Pkg1 = tilA*P0*tilA' + Q- Kt*Omegat*Kt'; % P(2|1)
    if min(eig(Pkg1))<0
        % this should not happen ideally (but does in practice)
        Pkg1 = Pkg1 + eye(n)*(0.1-min(eig(Pkg1)));
    end

    qlike = qlike +  log(det(Omegat)) + tres(1,:)*inv(Omegat)*tres(1,:)';
    Ts = size(ty,1);
    
    for t=2:Ts-1 % filter equations
        Omegat = tilC*Pkg1*tilC'+Omegas;
        iOm = inv(Omegat);
        tres(t,:)= ty(t,:) - xf'*tilC';
        Kt = (tilA*Pkg1*tilC'+tilK*Omegas)*iOm;
        xf = tilA*xf+ Kt*tres(t,:)';
        Pkg1 = (tilA * Pkg1 *tilA') +Q - (Kt*Omegat*Kt');
        Pkg1 = (Pkg1+Pkg1')/2;
        % update likelihood
        qlike = qlike + log(det(Omegat)) + tres(t,:)*inv(Omegat)*tres(t,:)';
    end
    
    % last year might be different.
    Omegat = tilC*Pkg1*tilC'+Omegas;
    
    % find out, how many obs are present?
    obs =sum(1-isnan(ty(Ts,:)));
    Omegat = Omegat(1:obs,1:obs);
    iOm = inv(Omegat);
    tres(Ts,:)=NaN;
    tres(Ts,1:obs)= ty(Ts,1:obs) - xf'*tilC(1:obs,:)';
    Kt = (tilA*Pkg1*tilC(1:obs,:)'+tilK*Omegas(:,1:obs))*iOm;
    xf = tilA*xf+ Kt*tres(t,1:obs)';
    Pkg1 = (tilA * Pkg1 *tilA') +Q - (Kt*Omegat*Kt');
    Pkg1 = (Pkg1+Pkg1')/2;
    % update likelihood
    qlike = qlike + log(det(Omegat)) + tres(t,1:obs)*inv(Omegat)*tres(t,1:obs)';
    
    
else
    P0 = 0*eye(n);
    % --- initialize the Kalman filter ---
    x0= zeros(n,1); % x(1|0)
    Abar = tilA-tilK*tilC;    
    % --- run the filter for the inverse  transfer function ---
    Ts = size(ty,1);
    xt = ltitr(Abar,tilK,ty(1:Ts-1,:),x0);
    tres = ty(1:(Ts-1),:)-xt*tilC';
    [tres,~] = est_initial_val(tres,Abar,tilK,tilC);
    
    xtt = Abar*xt(end,:)'+tilK*ty(Ts-1,:)';
    % last year might be different.    
    % find out, how many obs are present?
    obs =sum(1-isnan(ty(Ts,:)));
    
    tres(Ts,:)=NaN;
    tres(Ts,1:obs)= ty(Ts,1:obs) - xtt'*tilC(1:obs,:)';
    
    
    Omegat =  tres(1:(Ts-1),:)'*tres(1:(Ts-1),:)/(Ts-1);
    Omegat(1:obs,1:obs)= (Ts-1)/Ts*Omegat(1:obs,1:obs) + tres(Ts,1:obs)*tres(Ts,1:obs)'/Ts;
    % update likelihood
    qlike = qlike + Ts*log(det(Omegas)) + Ts*trace(Omegat*inv(Omegas));

end

% reshape residuals.
vr = inv(Es)*tres';
res = reshape(vr(:),s,Ts*S)';
res = res(1:T,:);
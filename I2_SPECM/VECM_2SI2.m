function [th_est,Q1,Q2,alpha,beta,zeta,eta,psi,Best,Gamma,Omega,et,beta_perp] = VECM_2SI2(y,p,rs,const)
% VECM_2SI2 implements the two step estimation procedure for I(2) VAR(p)
% processes od Paruolo (2000).
% 
% SYNTAX: [th_est,alpha,beta,zeta,eta,psi,B,Omega] = VECM_2SI2(y,p,s,m,rs);
%
% INPUTS:    y ... T x(s+m) matrix of observations
%            p ... integer; number of additional Delta^2 y_{t-j} terms.
%            rs ... 2x1 integer vector; no cols [beta_2,beta_1].
%            const ... indicator: use a constant or not?
%
% OUTPUT:   th_est ... theta structure
%           alpha,beta,eta, zeta, psi,B, Gamma ... matrices of I(2) VECM.
%           Omega ... sxs innovation variance estimate.
%
% AUTHOR: dbauer, 5.2.2020. 

% check inputs. 
[T,s]= size(y);

if (sum(rs)>s)
    error('VECM_2SI2: Something is wrong with dimensions rs.');
end

if p<0
    error('VECM_2SI2: lag length must be positive!');
end

% calculate first and second derivatives. 
dy = diff(y(:,1:s));
ddy = diff(dy);

Z0t = ddy; % Z0t = Delta^2 y_t
Z2t = y(2:end-1,:); % Z2t = y_{t-1} 
Z1t = dy(1:end-1,1:s); % Z1t = Delta y_{t-1}

% if constant -> concentrate out.
if const>0 
    Z0t = detrend(Z0t,0);
    Z1t = detrend(Z1t,0);
    Z2t = detrend(Z2t,0);
    ddy = detrend(ddy,0);
    Joh_j = 3;
else
    Joh_j = 5;
end

% build lags of ddy
if p>0
    Zst = ddy(1:end-1,:);
    for j = 1:(p-1) % add stationary terms
        Zst = [ddy(1:end-1,:),[NaN*ones(1,size(Zst,2));Zst(1:end-1,:)]];
    end
    Z0t = Z0t(p+1:end,:);
    Z1t = Z1t(p+1:end,:);
    Z2t = Z2t(p+1:end,:);
    Zst = Zst(p:end,:);
    
    % concentrate out stationary terms 
    R0t = Z0t - Zst*(Zst\Z0t);
    R1t = Z1t - Zst*(Zst\Z1t);
    R2t = Z2t - Zst*(Zst\Z2t);
else
    R0t = Z0t;
    R1t = Z1t;
    R2t = Z2t;
end


% Step 1: estimate alpha,beta from equation (3.4) without restriction of
% rank of Gamma.

% concentrate out R1t (Delta y_{t-1}).
S0t = R0t - R1t*(R1t\R0t);
S2t = R2t - R1t*(R1t\R1t);

[alpha,beta,Q1,~,~]= est_alpha_beta(S0t,S2t,rs(1),Joh_j);

[beta_q,beta_r] = qr(beta);
beta = beta_q(:,1:rs(1));
beta_perp = beta_q(:,rs(1)+1:end);
alpha = alpha*beta_r(1:rs(1),1:rs(1))';
[alpha_q,~]= qr(alpha);
alpha_perp = alpha_q(:,rs(1)+1:end);

% Step 2.1: estimate (3.9)
S0t = R0t*alpha_perp;
S1t = R1t*beta_perp;
Sbt = R1t*beta;

T0t = S0t - Sbt*(Sbt\S0t);
T1t = S1t - Sbt*(Sbt\S1t);

[zeta,eta,Q2,~]= est_alpha_beta(T0t,T1t,rs(2),Joh_j);

beta1 = beta_perp*eta;
[eta_q,~]= qr(eta);
beta2 = beta_perp*eta_q(:,rs(2)+1:end);

% Step 2.2: estimate psi in (3.8)
Ut = R0t*alpha*inv(alpha'*alpha)-R2t*beta;
Uta = Ut - S0t*(S0t\Ut);
U1t = R1t - S0t*(S0t\R1t);
psi = U1t\Uta;

% obtain estimator for delta from (2.5)
delta = psi'*beta2*inv(beta2'*beta2);

% last part: obtain estimate for B from (3.14)
Z0t = ddy; % Z0t = Delta^2 y_t
Z2t = y(2:end-1,:); % Z2t = y_{t-1} 
Z1t = dy(1:end-1,1:s); % Z1t = Delta y_{t-1}
if (p>1)
    Z0t = Z0t(p+1:end,:);
    Z1t = Z1t(p+1:end,:);
    Z2t = Z2t(p+1:end,:);
end
nt =[Z2t*beta+Z1t*beta2*delta',Z1t*[beta,beta1],Zst];

if const>0
    nt = [nt,ones(size(nt,1),1)];
    plus1=1;
else
    plus1=0;
end
rt = Z0t;

B = nt\rt;
et = rt - nt*B;
Omega = et'*et/(T-p-2);

% calculate VAR from VECM.
A = zeros(s,s*(p+2));
Best = B';
alpha= Best(:,1:rs(1));
Best(:,1:rs(1))=[];
xi1 = Best(:,[1:rs(1)]);
Best(:,1:rs(1))=[];
xi2 = Best(:,[1:rs(2)]);
Best(:,1:rs(2))=[];

Pi = alpha*beta';
Gamma = alpha*delta*beta2'+[xi1,xi2]*[beta,beta1]';


A(:,1:s)=2*eye(s)+Pi + Gamma;
A(:,(s+1):(2*s))=-eye(s)-Gamma; 

for j=3:(p+2)
    Psij = Best(:,(j-3)*s+[1:s]);
    A(:,(j-1)*s+[1:s]) = A(:,(j-1)*s+[1:s])+Psij;
    A(:,(j-2)*s+[1:s]) = A(:,(j-2)*s+[1:s])-2*Psij;
    A(:,(j-3)*s+[1:s]) = A(:,(j-3)*s+[1:s])+Psij;
end

th_est = theta_urs();
th_est.which = 'poly';
th_est.b = eye(s);
th_est.a = [eye(s),-A];
if plus1>0
    th_est.d = Best(:,end);
    th_est.m = plus1;
else
    th_est.d = zeros(s,1);
end
th_est.Omega = Omega;
th_est.ur = 'I(2)';

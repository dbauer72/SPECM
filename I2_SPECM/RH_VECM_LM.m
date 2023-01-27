function [th_init,Chat,Q1,Q2,alpha,beta,zeta,eta,psi,Best,Gamma] = RH_VECM_LM(y,Abar,B,rs,const);
% Ribarits_Hanzon proc for estimation of RH_VECM form for 2SI2 procedure.
%
% SYNTAX: [Chat,alpha,beta,zeta,eta,psi,B,Gamma,Omega] = RH_VECM_LM(y,Abar,B,rs,const);
%
% INPUT: y    ... Txs data matrix,
%        Abar ... nxn matrix A-BC.
%        B    ... nxs matrix
%        rs ... 2x1 integer vector; no cols [beta_2,beta_1].
%        const ... integer; indicator for constant term in VECM.
%
% OUTPUT: th_init ... theta stucture,
%         Chat    ... sxn optimal C.
%         alpha,beta,eta, zeta, psi,B, Gamma ... matrices of I(2) VECM.
%         Omega ... sxs innovation variance estimate.
%
%
% COMMENT: The restriction C*(I-Abar)^(-1)*B = I+alpha beta' is imposed.
%  Two step procedure in the light of the 2SI2 idea, but using
%  Lukas Matuschek idea of imposing zero restrictions based on initial
%  unrestricted estimates.
%
% AUTHOR: dbauer, 11.2.2020.

[T,s] = size(y);

if (sum(rs)>s)
    error('RH_VECM_LM: Something is wrong with dimensions rs.');
end

n = size(Abar,1);
ds(2) = rs(2);
ds(1) = s-sum(rs);


if (sum(ds)>n)
    error('RH_VECM_2SI2: Something is wrong with dimensions rs.');
end

[th_init,Chat,Q1,Q2,alpha,beta,zeta,eta] = RH_VECM_2SI2(y,Abar,B,rs,const,'n');


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

% --- filter ddy to obtain Z_{st} ----
Zst = zeros(T-2,n);
tB = B';
tA = Abar';
for t=2:size(ddy,1)
    Zst(t,:)=Zst(t-1,:)*tA + Z0t(t-1,:)*tB;
end;

Zst = Zst*inv(eye(n)-tA)^2*tA^2;

% Step 1: estimate alpha,beta from equation (3.4) without restriction of
% rank of Gamma.

% concentrate out Z1t (delta y_(t-1).
R0t = Z0t - Z1t*(Z1t\Z0t);
R2t = Z2t - Z1t*(Z1t\Z2t);
Rst = Zst - Z1t*(Z1t\Zst);

% set up restriction
[Q,R]=qr(beta);
beta_perp = Q(:,rs(1)+1:end);

R = inv(eye(n)-Abar)*B*beta_perp;
r = beta_perp;

% concentrate out Zst under the restriction on Chat.
tildeM = Rst'*Rst/T;
if min(eig(tildeM))<10^(-6) % regularize
    tildeM = tildeM + eye(size(Rst,2))*0.01;
end;
itM = inv(tildeM);
M0s = R0t'*Rst/T;
M2s = R2t'*Rst/T;

H11 = itM*(eye(n)-R*inv(R'*itM*R)*R'*itM);
H21 = beta_perp*inv(R'*itM*R)*R'*itM;

C0bar = M0s*H11 + H21;
C1bar = M2s*H11; % - H21;

% generate S_{0t} (R_{0t} corrected for estimated contribution of
% stationary vars)
S0t = R0t - Rst*C0bar';

% generate R_{1t} (R_{1t} corrected for estimated contribution of
% stationary vars)
S2t = R2t - Rst*C1bar';
% --- rank restricted regression for alpha and beta ---
[alpha,~,Q1]= est_alpha_beta(S0t,S2t,rs(1),5); % non constant or linear trend.

Ch = C0bar - alpha*beta'*C1bar; 
Pi = eye(s) - Ch*inv(eye(n)-Abar)*B;
alpha = Pi*beta*inv(beta'*beta); 

% refactor
%[beta_q,beta_r] = qr(beta);
%beta = beta_q(:,1:rs(1));
%beta_perp = beta_q(:,rs(1)+1:end);
%alpha = alpha*beta_r(1:rs(1),1:rs(1))';

% estimate of Calpha
Calpha = inv(alpha'*alpha)*alpha'*Ch;

% calculate alpha_perp
[alpha_q,~]= qr(alpha);
alpha_perp = alpha_q(:,rs(1)+1:end);

% Step 2.1: estimate (3.9) of Paruolo, starting from unadjusted
% vars
R0t = Z0t*alpha_perp;
R1t = Z1t*beta_perp;
Rbt = Z1t*beta;

S0t = R0t - Rbt*(Rbt\R0t);
S1t = R1t - Rbt*(Rbt\R1t);
Sst = Zst - Rbt*(Rbt\Zst);

% concentrate out hat C_perp tilde x_t
% set up constraints
% calculate eta_perp
[eta_q,~]=qr(eta);
eta_perp = eta_q(:,(rs(2)+1):end);
Ra = [inv(eye(n)-Abar)*B*beta_perp,inv(eye(n)-Abar)^2*B*beta_perp*eta_perp];
ra = [alpha_perp'*beta_perp,zeros(s-rs(1),s-rs(1)-rs(2))];

% --- linearly restricted regression ---
tildeM = Sst'*Sst/T;
if min(eig(tildeM))<10^(-6) % regularize
    tildeM = tildeM + eye(size(Rst,2))*0.01;
end;
itM = inv(tildeM);
M0s = S0t'*Sst/T;
M1s = S1t'*Sst/T;

RiMR = Ra'*itM*Ra;
if min(abs(eig(RiMR)))<0.00001
    RiMR = RiMR + eye(size(RiMR,1))*0.00001;
end;

H11 = itM*(eye(n)-Ra*inv(RiMR)*Ra'*itM);
H21 = inv(RiMR)*Ra'*itM;

C0abar = M0s*H11 + ra*H21;
C1abar = M1s*H11; % + H21((s+1):end,:);

% generate S_{0t} (R_{0t} corrected for estimated contribution of
% stationary vars)
T0t = S0t - Sst*C0abar';

% generate R_{1t} (R_{1t} corrected for estimated contribution of
% stationary vars)
T1t = S1t - Sst*C1abar';

% --- rank restricted regression for zeta and eta ---
[zeta,~,Q2,~]= est_alpha_beta(T0t,T1t,rs(2),5);

% --- complete estimation of Chat ---
Cperp = C0abar - zeta*eta'*C1abar;

Chat = [alpha,alpha_perp]*[Calpha;Cperp];

% --- estimate remaining quantities ---
beta1 = beta_perp*eta;
[eta_q,~]= qr(eta);
beta2 = beta_perp*eta_q(:,rs(2)+1:end);

% Step 2.2: estimate psi in (3.8)
Ut = (Z0t-Zst*Chat')*alpha*inv(alpha'*alpha)-Z2t*beta;
Sut = (Z0t-Zst*Chat')*alpha_perp;
Uta = Ut - Sut*(Sut\Ut);
U1t = Z1t - Sut*(Sut\Z1t);
psi = U1t\Uta;

% obtain estimator for delta from (2.5)
delta = psi'*beta2*inv(beta2'*beta2);

% last part: obtain estimate for B from (3.14)
nt =[Z2t*beta+Z1t*beta2*delta',Z1t*[beta,beta1]];
rt = Z0t-Zst*Chat';

Bb = nt\rt;
et = rt - nt*Bb;
Omega = et'*et/size(et,1);

Best = Bb';
%alpha= Best(:,1:rs(1));
Best(:,1:rs(1))=[];
xi1 = Best(:,[1:rs(1)]);
Best(:,1:rs(1))=[];
xi2 = Best(:,[1:rs(2)]);
Best(:,1:rs(2))=[];

Pi = alpha*beta';
Gamma = alpha*delta*beta2'+[xi1,xi2]*[beta,beta1]';

ds = [s-sum(rs),rs(2)];
%[Ae,Be,Ce] = convert_canform_I2(Abar+B*Chat,B,Chat,rs);
Ae = Abar+B*Chat;
Be = B;
Ce = Chat;

th_init = theta_urs();
th_init.A = Ae;
th_init.K = Be;
th_init.C = Ce;
th_init.which = 'SS';
th_init.urs = [ds];
th_init.Omega = eye(s);
th_init.ur = 'I(2)';

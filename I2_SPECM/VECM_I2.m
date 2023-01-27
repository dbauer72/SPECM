function [th_est,alpha_sw,beta_sw,Pi_sw,Best,Gamma_sw,Omega,res,Q1,Q2] = VECM_I2(y,p,rs,const);
% VECM form estimation of AR systems in the I(2) case following Johansen's
% two step estimator in the VECM formulation.
%
% SYNTAX: [LL,res,th] = VECM_I2(y,p,rs,Joh_j);
%
% INPUT: y    ... Txs data matrix,
%        p    ... lag order of VECM Delta terms.
%        rs    ... 2x1 int; coint-rank for the two steps.
%        Joh_j ... deterministics of Johansen. 
%
% OUTPUT: LL   ... minimized log likelihood.
%         res  ... Txs residuals
%         th   ... theta structure containing the ARX system estimated. 
%         u    ... Txm real matrix of deterministic terms included.
%
% REMARK: Joh_j is an integer containing the specification of the
% characteristics according to Johansen (1997). 
%  Joh_j:    1 ... no restriction. 
%
% AUTHOR: dbauer, 21.8.2020 (based on procedure from Yuanyuan Li). 

% prepare regressions 
if nargin< 4
    Joh_j = 5; % no deterministic terms. 
else
    if const>0
        Joh_j = 3;
    else
        Joh_j = 5;
    end;
end;
    
[T,s] = size(y);

if p<0 % negative lag length implies AIC estimation of order.
    p = aicest(y,s,abs(p));
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
end

% build lags of ddy
if p>0
    Zst = ddy(1:end-1,:);
    for j = 1:(p-1) % add stationary terms
        Zst = [ddy(1:end-1,:),[NaN*ones(1,size(Zst,2));Zst(1:end-1,:)]];
    end;
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


% first step RRR: d2Y = Pi dLY + stationary terms.
% concentrate out R2. 
R0pt = R0t - R1t*(R1t\R0t);
R2pt = R2t - R1t*(R1t\R2t);

[alphahat,betahat,Q1]= est_alpha_beta(R0pt,R2pt,rs(1),Joh_j);

% 2nd step: RRR of alpha_perp' R0t on R1t_perp = R2t beta_perp
[beta_q,beta_r] = qr(betahat);
beta_h = beta_q(:,1:rs(1));
beta_o_h = beta_q(:,rs(1)+1:end);
alpha_h = alphahat*beta_r(1:rs(1),1:rs(1))';
[alpha_q,~]= qr(alpha_h);
alpha_o_h = alpha_q(:,rs(1)+1:end);

% Step 2.1: estimate (3.9)
S0t = R0t*alpha_o_h;
S1t = R1t*beta_o_h;
Sbt = R1t*beta_h;

T0t = S0t - Sbt*(Sbt\S0t);
T1t = S1t - Sbt*(Sbt\S1t);

[zeta_h,eta_h,Q2]= est_alpha_beta(T0t,T1t,rs(2),Joh_j);

[eta_q,eta_r] = qr(eta_h);
eta_h = eta_q(:,1:rs(2));
eta_o_h = eta_q(:,rs(2)+1:end);
%zeta_h = zeta_h*eta_r(1:rs(1),1:rs(1))';
%[zeta_q,~]= qr(zeta_h);
%zeta_o_h = zeta_q(:,rs(1)+1:end);

% after initialization, start switching. 

beta1_h = beta_o_h * inv(beta_o_h'*beta_o_h) * eta_h;
beta2_h = beta_o_h * eta_o_h;
        
% initialization of tau
tau = [beta_h beta1_h];
tau_o = beta2_h;

[nrtau,nctau]=size(tau); 

%%%%%%%%%%%%%%%%%
%%% switching %%%
%%%%%%%%%%%%%%%%%


di = 1;  % norm(Pi_h-Pi_Tr,1);
Pi_sw = zeros(s);
Gamma_sw = zeros(s);
  
r_sw = rs(1);
s_sw = rs(2);

it = 1; 
while (di > 0.000001)&&(it<20)
        
     % Given tau and tau_o, estimate other parameters.
     % RRR d^2 X_t on (tau'X_t-1; tau_o'dX_t-1) corrected for tau'dX_t-1
        R0_sw = R0t;
        R2_sw = [ R2t*tau , R1t*tau_o];
        R1_sw = R1t*tau;
        
        R0st = R0_sw - R1_sw*(R1_sw\R0_sw);
        R2st = R2_sw - R1_sw*(R1_sw\R2_sw);

        [alpha_sw,beta_sw,~]= est_alpha_beta(R0st,R2st,rs(1),Joh_j);
        [alpha_q,~]= qr(alpha_sw);
        alpha_o_sw = alpha_q(:,rs(1)+1:end);

        res_sw = R0st - R2st*beta_sw*alpha_sw';         
        Omega_sw = (res_sw'*res_sw)/(T-p);
        
        % rho, alpha_o_Omega
        rho_sw = beta_sw(1:(r_sw+s_sw),:);  % rho = [I_r; 0]: (r+s) x r.
        alpha_o_Omega_sw = Omega_sw*alpha_o_sw*inv(alpha_o_sw'*Omega_sw*alpha_o_sw);      
%         delta_sw = beta_sw(r_sw+s_sw+1:end,:)'; 
        
        % kappa 
        y_z = R0_sw*alpha_o_sw;          
        kappa_sw = (R1_sw\y_z);        % kappa_sw: kappa', (p-r) x (r+s). 
        
        Omega_1_hi = (alpha_h'*inv(Omega_sw)*alpha_h);
        Omega_2_h = alpha_o_h'*Omega_sw*alpha_o_h;
        
        % equation on p. 451 of Joha (1997; Scandinavian journal of
        % Statistics). 
        
        % left hand side 
        aod2 = [R0t*alpha_o_h,Z1t];
        Xt1_perp = Z2t - aod2*(aod2\Z2t);
        
        E = rho_sw*(Omega_1_hi)*alpha_h'*Z0t'*Xt1_perp + kappa_sw*inv(Omega_2_h)*alpha_o_h'*Z0t'*Z2t; 
        
        %right hand side
        A = rho_sw*(Omega_1_hi)*rho_sw';
        B = Xt1_perp'*Xt1_perp;
        C = kappa_sw*inv(Omega_2_h)*kappa_sw';
        D = Z2t'*Z2t;
        
        % solve equation 
        vec_tau = inv(kron(B',A)+kron(D,C'))*E(:);

        tau_temp = reshape(vec_tau, nctau,nrtau)';  % tau, p x (r+s)
        [tau_q,~]= qr(tau_temp);
        
        tau_temp = tau_q(:,1:nctau);
        % psi
        yy = R0t*inv(Omega_sw)*alpha_h*Omega_1_hi - Z1t*tau_temp*rho_sw;
        ZZ = R1t; 
        psi_sw = ZZ\yy;
        
        % Pi, Gamma
        Pi_sw_temp = alpha_sw * rho_sw' * tau_temp';
        Gamma_sw_temp =   alpha_sw*psi_sw' + alpha_o_Omega_sw*kappa_sw'*tau_temp';
        
        di = max(norm( Pi_sw_temp - Pi_sw,1), norm(Gamma_sw_temp - Gamma_sw,1)) ;
        
        Pi_sw = Pi_sw_temp;
        Gamma_sw = Gamma_sw_temp;
   
        tau = tau_temp;
        
        % update of tau_o
        [Q,~]=qr(tau);
        tau_o = Q(:,(r_sw+s_sw)+1:end); 
            
        % increase count to eliminate infinite loops. 
        it = it+1; 
end;


%%%%%%%%%%%%%%%%%%%%
%%%% end switch %%%%
%%%%%%%%%%%%%%%%%%%%

% estimate the remaining parameters
yy = Z0t - Z2t*Pi_sw' - Z1t * Gamma_sw'; 

Best = (Zst\yy)';

res = yy - Zst*Best'; 
Omega = res'*res/(T-p-2); 

if Joh_j==3
    Z0t = ddy; % Z0t = Delta^2 y_t
    Z2t = y(2:end-1,:); % Z2t = y_{t-1} 
    Z1t = dy(1:end-1,1:s); % Z1t = Delta y_{t-1}    
    Z0t = Z0t(p+1:end,:);
    Z1t = Z1t(p+1:end,:);
    Z2t = Z2t(p+1:end,:);
    
    res2t = Z0t - Z2t*Pi_sw' - Z1t * Gamma_sw' - Zst*Best';
    D_h = mean(res2t)';
else
    D_h = zeros(s,1);
end

% convert VECM to VAR formulation.
A = zeros(s,s*(p+2));
A(:,1:s)=2*eye(s)+Pi_sw + Gamma_sw;
A(:,(s+1):(2*s))=-eye(s)-Gamma_sw; 

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
th_est.d = D_h;

th_est.Omega = Omega;




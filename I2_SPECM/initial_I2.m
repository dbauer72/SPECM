function th_init = initial_I2(y,rs,n,f);

% implements initialization procedure of Li, Y. and Bauer Econometrics
% (2020). 

[T,p]= size(y); 
% Step 1: estimate 2SI2 VAR using f differenced stationary terms.
[th_est,Q1,Q2,alpha,beta,zeta,eta,psi,Best,Gamma,Omega,et,beta_perp] = VECM_2SI2(y,f,rs,0);
th_est.m = 1;

% obtain impulse response sequence. 
th_ss = poly2ss_th(th_est);
Ae = th_ss.A;
Ke = th_ss.K;
Ce = th_ss.C; 

Ait = eye(size(Ae,1));

for j=1:(2*f+1)
    Pij(:,:,j)=-Ce*Ait*Ke;
    Ait = Ait*(Ae-Ke*Ce);
end

% Step 2: obtain (Abar,B,D) from Pij. 
H = zeros(p*f,p*f);
for i=1:f
    for j=1:f
        H((i-1)*p+[1:p],(j-1)*p+[1:p])=Pij(:,:,(i+j-1));
    end
end

[U,S,V]= svd(H);
Of = U(:,1:n); 
B = S(1:n,1:n)*V(1:p,1:n)';
Abar = Of(1:(end-p),:)\Of(p+1:end,:);

% Step 3: project row of eta'
Beta = B*beta_perp;

if size(Beta,1)<size(Beta,2)
    bbb = Beta';
    etatilde = ( bbb*(bbb\eta) );
else
    etatilde = eta;
end;


% Step 4: Solve equations 
R2 = inv(eye(n)-Abar)^2*Beta;
r2 = -zeta*etatilde'; 

Cdag0 = (R2*inv(R2'*R2)*r2')';

% orthogonal direction 
[Q,R] = qr(R2);
R_perp = Q(:,p-rs(1)+1:end);

prm = size(Cdag0,1);
rp = size(R_perp,2);

if size(R_perp,2)>0

    options = optimoptions('fminunc','display','iter','MaxIterations',1000);
    options.MaxFunctionEvaluations = 20000;

    [cdel,fval,exitflag] = fminunc(@(x) cal_crit_initial(x,R_perp,Cdag0,Best,Abar,B,beta),randn(prm*rp,1),options);

    cdel = randn(prm*rp,1);

    [cr,Cest] = cal_crit_initial(cdel,R_perp,Cdag0,Best,Abar,B,beta);    
else 
    alpha_perp = Cdag0*inv(eye(n)-Abar)*B;
    Cest = inv(alpha_perp)*Cdag0;
end;

% Step 5: transform to canonical form. 
ds = [p-sum(rs),rs(2)];
[Ae,Be,Ce] = convert_canform_I2(Abar+B*Cest,B,Cest,rs); 
th_init = theta_urs();
th_init.A = Ae;
th_init.K = Be;
th_init.C = Ce;
th_init.which = 'SS';
th_init.urs = [ds];
th_init.Omega = eye(p);
th_init.ur = 'I(2)';
 



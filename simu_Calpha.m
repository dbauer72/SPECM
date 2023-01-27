function [alphahat,alphahat_1,alphahat_2,alphahat_3,alphahat_v,alpha] = simu_Calpha(A,K,C,beta,T,M);
% simulate and estimate using Calpha method. 

% initialisation 
Joh_j = 5; 
%T = 2000;
[s,n]= size(C);
restrict.det_res = 0;

r = size(beta,2);

% calulcate alpha 
Abar = A-K*C;
Pi = eye(s) - C*inv(eye(n)-Abar)*K;
alpha = -Pi*beta*inv(beta'*beta);
baralpha = alpha*inv(alpha'*alpha);
h = waitbar(0,'Please wait ...');
for m= 1:M
    waitbar(m/M,h);
    % generate y 
    x = zeros(n,1);
    u = randn(T,n);
    y = u*0; 

    for j=1:T
        y(j,:) = x'*C'+ u(j,:);
        x = A*x + K*u(j,:)';
    end;

    % estimate system 
    [k,sig,kbc,sigbc,phi,sigphi,thar] = aicest(y,s,max(s,floor(sqrt(T))));
    [th,Ah,Kh,Ch,Omega] = CCA(y,n,2*k,2*k,0);
    Abar = Ah-Kh*Ch;
    [~,~,betahat,~] = RH_vecm(y,Abar,Kh,size(beta,2),'n',Joh_j);

    % estimate alpha: without correction
    correct = 0;
    [alpha_m,Va] = Calpha(y,Abar,Kh,betahat,Joh_j,correct);
    % align with true to avoud problems with normalisation.
    [alpha_m,Va] = standardize_alpha(alpha_m,baralpha,Va);
    % test first row = 0.
    tval = alpha_m(1,:)*inv(Va(1:r,1:r))*alpha_m(1,:)';
    alphahat(m,:)= [(alpha_m(:))',tval];

    % estimate alpha: with correction
    correct = 1;
    [alpha_m,Va] = Calpha(y,Abar,Kh,betahat,Joh_j,correct);
    % align with true to avoud problems with normalisation.
    [alpha_m,Va] = standardize_alpha(alpha_m,baralpha,Va);
    tval = alpha_m(1,:)*inv(Va(1:r,1:r))*alpha_m(1,:)';
    alphahat_1(m,:)= [(alpha_m(:))',tval];

    % first optimize stationary state space model
    [result,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM_I1(y,s,n,0,n,0,restrict,th);
    [~,~,betahat,~] = RH_vecm(y,Ac-Kc*Cc,Kc,size(beta,2),'n',Joh_j);

    [alpha_m,Va] = Calpha(y,Ac-Kc*Cc,Kc,betahat,Joh_j,0);
    % align with true to avoud problems with normalisation.
    [alpha_m,Va] = standardize_alpha(alpha_m,baralpha,Va);
    tval = alpha_m(1,:)*inv(Va(1:r,1:r))*alpha_m(1,:)';
    alphahat_2(m,:)= [(alpha_m(:))',tval];

    [alpha_m,Va] = Calpha(y,Ac-Kc*Cc,Kc,betahat,Joh_j,1);
    % align with true to avoud problems with normalisation.
    [alpha_m,Va] = standardize_alpha(alpha_m,baralpha,Va);
    tval = alpha_m(1,:)*inv(Va(1:r,1:r))*alpha_m(1,:)';
    alphahat_3(m,:)= [(alpha_m(:))',tval];   

    % use VECM approach.
    [~,alpha_m] = VECM_I1(y,k,size(beta,2),Joh_j);
    [alpha_m,Va] = standardize_alpha(alpha_m,baralpha,Va);
    tval = alpha_m(1,:)*inv(Va(1:r,1:r))*alpha_m(1,:)';
    alphahat_v(m,:)= [(alpha_m(:))',tval];
end

close(h); 

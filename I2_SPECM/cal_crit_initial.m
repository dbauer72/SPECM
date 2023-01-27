function [crit,Cest] = cal_crit_initial(cdel,R_perp,Cdag0,Pij,Abar,B,beta);
% cal_crit_initial calculates the least squares distance of Pij to the
% corresponding estimate from (A,B,C).
%
% SYNTAX: [crit,Cest] = cal_crit_initial(cdel,Cdag0,Pij,Abar,B,beta);
%
% INPUT: cdel ... vector of parameters.
%        Cdag0 ... p-r x n real matrix; particular solution to equation(3).
%        Pij   ... p x p*f real matrix of estimates of impulse response.
%        Abar  ... n x n matrix of inverse system.
%        B     ... n x p Kalman gain of inverse system.
%        beta  ... p x r real matrix of cointegration space.
%
% OUTPUT: crit ... real; norm of distance.
%         Cest ... corresponding estimated C matrix.
%
% remark: implements Step 4 of the initialization procedure of Li and Bauer (2020).
%
% AUTHOR: dbauer, 1.9.2020. 

% solve second equation 
[pmr,p]=size(Cdag0);
[n,s]=size(B);
r = size(beta,2); 

rp = size(R_perp,2);

[~,p]=size(B);
r = p-pmr;
Cdel = reshape(cdel,pmr,rp);

alpha_perp = (Cdag0+Cdel*R_perp')*inv(eye(n)-Abar)*B;

% set up least squares criterion
Cf = zeros(n,size(Pij,2));

f = round(size(Pij,2)/p);
Cf(:,1:p)=B;
for j=2:f
    Cf(:,(j-1)*p+[1:p])= Abar*Cf(:,(j-2)*p+[1:p]);
end;

Cf = inv(eye(n)-Abar)^2*Cf;

    
if r>0 % are there stationary components?
    [Q,R]= qr(alpha_perp');
    alpha_o = Q(:,(end-r+1):end);

    % set up first equation 
    R1 = [-alpha_o';inv(eye(n)-Abar)*B];
    r1 = beta'; 

    [Q,R]= qr(R1);
    R1_perp = Q(:,p+1:end);

    TC = (R1*inv(R1'*R1)*r1')';
    Ca = TC(:,r+1:end);
    RCadel = R1_perp(r+1:end,:)';
    alphatil = alpha_o*TC(1:r,1:r)';

    

    % calculate difference.
    lhs = alphatil'*Pij - Ca*(Abar^2)*Cf;
    rhs = RCadel*Cf;

    Cadel = lhs*rhs'*inv(rhs*rhs');

    % estimate
    TCo = TC(1:r,1:r)+Cadel*R1_perp(1:r,:)';
    alphatil = alpha_o*TCo';
    
    Ta = [alphatil';alpha_perp];
    Cest = inv(Ta)*[Ca+Cadel*RCadel;Cdag0+Cdel*R_perp'];

else
    Cest = Cdag0+Cdel*R_perp';
end;


% criterion function. 
crit = norm(Pij - Cest*(Abar^2)*Cf);





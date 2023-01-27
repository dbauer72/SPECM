function pred = cal_prediction(y,dt,th,H);
% calculate the prediction for the time series y using the model th,
% where the deterministics are given by dt.
% Predictions are calculated at the end of the sample for up to H time
% periods ahead. 
%
%  SYNTAX: pred = cal_prediction(y,dt,th,H);
% 
% INPUT: y ... Txs time series
%        dt ... Txm deterministics
%        th ... theta structure.
%        H ... time horizon.
%
% OUTPUT: pred ... Txs predictions (out of sample; filled up with one step
% preds before).
%
% REMARKS: uses the Kalman filter beforehand. 
% AUTHOR: dbauer, 13.12.2019.

A = th.A;
K = th.K;
C = th.C;
D = th.D; 
B = th.B; 
Omega = th.Omega;
Q = K*Omega*K'; 

tilA = A;
tilK= K;
tilC=C;
    
m = size(dt,2);
[T,s] = size(y);

ty = y - dt*D';
n = size(A,1);

% --- initialize the Kalman filter ---
x0= zeros(n,1);P0 = zeros(n,n);
tres = ty*0;
tres(1,:)=ty(1,:);
Omegat= tilC*P0*tilC'+Omega;
Kt = (tilA*P0*tilC'+tilK*Omegat)*inv(Omegat);
xf = tilA*x0 + Kt*tres(1,:)';
Pkg1 = tilA*P0*tilA' + Q- Kt*Omegat*Kt'; % P(2|1)

for t=2:T-H % filter equations unitl T-H. 
    Omegat = tilC*Pkg1*tilC'+Omega;
    iOm = inv(Omegat);
    tres(t,:)= y(t,:) - xf'*tilC' - dt(t,:)*D';
    Kt = (tilA*Pkg1*tilC'+tilK*Omega)*iOm;
    xf = tilA*xf+ Kt*tres(t,:)' + B*dt(t,:)';
    Pkg1 = (tilA * Pkg1 *tilA') +Q - (Kt*Omegat*Kt');    
    Pkg1 = (Pkg1+Pkg1')/2; 
end

pred = y-tres;  
for j=1:H
    pred(T-H+j,:)= (tilC*xf + D*dt(T-H+j,:)')';
    xf = tilA*xf + B*dt(T-H+j,:)';
end


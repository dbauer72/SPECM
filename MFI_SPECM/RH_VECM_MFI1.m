function [LL,th_est,Chat] = RH_VECM_MFI1(y,S,Abar,B,estimate,Joh_j);
% Ribarits_hanzon idea for estimation of RH_VECM_MFI1 form in the MFI1
% case.
%
% SYNTAX: [LL,th_est,Chat] = RH_VECM(y,Abar,B,r,rest,Joh_j);
%
% INPUT: y    ... Txs data matrix,
%        S    ... integer; number of seasons;
%        Abar ... nxn matrix A-BC.
%        B    ... nxs matrix
%        estimate ... structure of cell arrays; contains initial estimates
%                of alphas, betas. Alternatively unit root structure.  
%        Joh_j ... integer; indicating the version of Johansen's
%        deterministics. 
%
% OUTPUT: LL   ... minimized log likelihood.
%         th_est ... theta structure of state space system.
%         Chat ... sxn optimal C.
%
% REMARK: the procedure implements the 'restricted' estimate holding the
%        betas fixed as described in the dissertation of Lukas Matuschek.
%
% AUTHOR: dbauer, 9.10.2020. 

[T,s]=size(y); 
% prepare regressions 
if nargin< 5
    Joh_j = 5; % no deterministic terms. 
end;

if nargin<4 
    error('RH_VECM_MFI: No urs supplied!');
end;

if ~isa(estimate,'vecm_mfi') % estimate does not yet contain the estimated Pi matrices. 
    urs = estimate;
    
    % select k with AIC
    k = aicest(y,s,floor(sqrt(T)));
    if k<S
        k=S;
    end;
    % perform AR approximation and estimation. 
    [~,~,~,~,estimate] = VECM_MFI1(y,s,S,urs,Joh_j);    
else
    urs = estimate.urs;    
end

n = size(Abar,1);
% steps: calculate tilde x_t: filtered version of Delta y.

[pX,pXk,Xm2,~,urs_full] = cal_reg_urs(y,S,s); 

% deal with deterministics 
[pXk,Xm2] = cal_reg_det(S,pXk,Xm2,urs_full,Joh_j);

pXk = zeros(size(pXk,1),0); % for the state space case, different stat. regressors are needed. 

% find out, which unit roots are present. 
urs(:,4)=0;
for j=1:size(urs_full,1)
    abu = abs(urs(:,1)-urs_full(j,1));
    [mi,mind]=min(abu);
    if mi>0.0001 %unit root not contained
        XX = squeeze(Xm2(:,:,j));
        if urs_full(j,2)==1
            pXk = [pXk,XX];
        else
            pXk = [pXk,XX(:,1:s)];
            if Joh_j == 4 % add regressor to coint. relation.
                pXk = [pXk,XX(:,2*s+1)]; 
            end;
            if (Joh_j<4)&&(urs_full(j,1)>0.0001)
                pXk = [pXk,XX(:,2*s+1)]; 
            end;
        end
    else % unit root contained -> remember index. 
        urs(mind,4)=j;
    end
end


in = find(urs(:,4)>0); % only these roots are contained. 
% compile the matrix of regressors for these: 

lin = length(in);
Xtk = zeros(T,0);

ind_k = cell(lin,1); 

cur = 1;
J = zeros(s,0); 
BkJ = zeros(n,0); 
G = zeros(0,0);

for j=1:lin
    jj = in(j);
    XX = squeeze(Xm2(:,:,jj));
    if urs_full(jj,2)==1 % complex unit root
        zk = exp(pi*sqrt(-1)*urs_full(jj,1));
        uBk = zk*inv(eye(n)-zk*Abar)*B;
        bbeta = estimate.beta{j};
        %bbeta = [real(beta),-imag(beta);imag(beta),real(beta)];
        rk = size(bbeta,2);
        curJ = size(J,2);
        Xtk = [Xtk,XX*bbeta];        
        G(end+[1:rk],1:curJ+2*s)= 0; 
        J = [J,eye(s),zeros(s,s)];
        BkJ = [BkJ,real(uBk),imag(uBk)];
        G(end-rk+[1:rk],curJ+1:end)=-bbeta';
    else
        zk = exp(pi*sqrt(-1)*urs_full(jj,1));
        uBk = zk*inv(eye(n)-zk*Abar)*B;
        bbeta = estimate.beta{j};
        Xtk = [Xtk,XX(:,1:s)*bbeta];        
        rk = size(bbeta,2);
        curJ = size(J,2);
        G(end+[1:rk],1:curJ+s)= 0; 
        J = [J,eye(s)];
        BkJ = [BkJ,real(uBk)];
        G(end-rk+[1:rk],curJ+1:end)=-bbeta';
        
        %if Joh_j == 4 % add regressor to coint. relation.
        %    Xtk = [Xtk,XX(:,2*s+1)];
        %end;
    end
    %if (Joh_j<4)&&(urs_full(jj,1)>0.0001)
    %   Xtk = [Xtk,XX(:,2*s+1)];
    %end;   
    
    ind_k{j} = [cur,size(Xtk,2)];
    cur = size(Xtk,2)+1;
end

% compose G and Bk
G = [BkJ;G]; 

% no rename to conform to the notation of Lukas Matuschek.
Z0t = pX;
vt = zeros(T-S,n);

tB = ((-inv(eye(n)-Abar^S)*Abar^S)*B)';
tA = Abar';
vt(1,:)=Z0t(1,:)*tB;

for t=2:T-S
    vt(t,:)=vt(t-1,:)*tA + Z0t(t-1,:)*tB;
end;


Z2t = vt;
%Z0t = Z0t(S+1:end,:);
%Xtk = Xtk(S+1:end,:); 

% pXk contains the deterministic terms and filtered terms not corresponding
% to unit roots. 
ZpK = pXk; %(S+1:end,:);

% calculate residuals with respect to Zpk;
R2t = Z2t - ZpK*(ZpK\Z2t);
R0t = Z0t - ZpK*(ZpK\Z0t);
Rtk = Xtk - ZpK*(ZpK\Xtk);

% all remaining regressors 
Vt = [R2t,Rtk]; 

% now estimate using the restriction. 

MM = [ Vt'*Vt,G;G',zeros(size(G,2),size(G,2))];
HH = inv(MM);

Chat = [ R0t'*Vt, J ]*HH(:,1:n);

Hhat = [ R0t'*Vt, J ]*HH(:,1:size(Vt,2)); 
resi = R0t - Vt*Hhat';
Omegai = resi'*resi/T; 


Ahat = Abar + B*Chat;

% write into structure
th_est.which = 'SS';
th_est.A = Ahat;
th_est.K= B;
th_est.C = Chat;
th_est.D = zeros(s,0); 
th_est.Omega = Omegai;

th_est.ur = 'MFI(1)';
th_est.urs = urs;

th_est = trans_can_form_MFI(th_est,urs);

[LL,res] =  cal_quasi_like_theta_MFI1(th_est,y,s,S,n,urs,-1);
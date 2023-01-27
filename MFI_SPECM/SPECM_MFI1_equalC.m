function  [result,thc,Ac,Kc,Cc,Omegac,thi_approx] = SPECM_MFI1_equalC(z,s,S,n,urs,nmax,Pbull,restrict,thi);
% SPECM_MFI1_equalC optimizes the Gaussian pseudo likelihood starting from an inital
% estimated system.
% 
% SYNTAX: [result,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM_MFI1_equalC(z,s,S,n,urs,nmax,Pbull,restrict,thi);
%
% INPUT:  z ... Tx(s+m) data: [y_t,d_t]
%         s ... integer; output dimension
%         S ... integer; number of seasons (quarterly data: S=4).
%         n ... state order to use
%         urs ... MFI1 unit root structure; number of common trends for
%                   each unit root. 
%         nmax ... maximum for AR order selection.
%         Pbull ... indicator; if Pbull>0, state is started in stationary
%                   distribution.
%         restrict ... structure specifying restrictions.
%         thi    ... use initial estimate.
%
% OUTPUT: result ... structure containing the estimation results.
%         th     ... theta estimate
%         (Ac,Kc,Cc) ... state space system in canonical form.
%         Omegac ... noise variance estimate
%         lle   ... minimizing value of the unrestricted quasi scaled log
%                  likelihood (over M_n). 
%
% AUTHOR: dbauer, 15.12.2021

% make sure equal is activated in restrict.
restrict.equal= 1;
% convert to canonical form
thi = trans_can_form_MFI(thi,urs);
A=  thi.A;
K = thi.K;
C = thi.C;
ck = urs(1,3); 

% approximate Cks 
Ck = approx_equalC(C,urs); 
cur=0;
for j=1:size(urs,1)
    if urs(j,2)==0 % real root 
        Ckq = Ck'*C(:,cur+[1:ck]);
        C(:,cur+[1:ck])=Ck;
        K(cur+[1:ck],:) = Ckq*K(cur+[1:ck],:);
        cur = cur+ck;
    else
        Ckq = Ck'*(C(:,cur+[1:ck])+sqrt(-1)*C(:,cur+ck+[1:ck]));
        [qc,r]=qr(Ckq);
        C(:,cur+[1:2*ck]) = [real(Ck*qc),imag(Ck*qc)];
        K(cur+[1:ck],:) = r*K(cur+[1:ck],:);
        K(cur+ck+[1:ck],:) = r*K(cur+ck+[1:ck],:);
        cur = cur+2*ck;
    end
end

thi_approx= thi;
thi_approx.K = K;
thi_approx.C = C;

% calculate parameter vector
parami = th2param_MFI1(thi_approx,urs,0,restrict);
restrict.scale = ones(length(parami),1);

parc = est_cal_like_hess_MFI1(z,s,S,n,urs,thi_approx,restrict);
param = parc;

% random improvement of I(1) systems
result = compile_results_MFI1(parc,s,S,n,urs,restrict,z(:,1:s),z(:,s+1:end),Pbull,nmax);
llc_c = result.deviance;
llc=llc_c+1;
it = 0;
while ((it<5)&&(llc_c<llc))
    if it>0
        disp(sprintf('Running random improve %2d. Old value: %4.2f -> new value: %4.2f.',it,llc,llc_c));
    end
    llc = llc_c;
    result = random_improve_call(result,0.001,5);
    llc_c =result.deviance;
    it = it+1;
end;

% write out results.
parc = result.param;
result = compile_results_MFI1(parc,s,S,n,urs,restrict,z(:,1:s),z(:,s+1:end),Pbull,nmax);
result.theta.ur='MFI(1)';
thc = result.theta;
m = size(z(:,s+1:end),2);
thc.B = zeros(n,m);
Ac = thc.A;
Kc = thc.K;
Cc = thc.C;
Omegac = thc.Omega;

function  [result,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM_MFI1(z,s,S,n,urs,nmax,Pbull,restrict,thi);
% SPECM optimizes the Gaussian pseudo likelihood starting from an inital
% CVA estimated system.
% 
% SYNTAX: [result,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM_MFI1(z,s,S,n,urs,nmax,Pbull,restrict);
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
%         thi    ... use initial estimate rather than CCA.
%
% OUTPUT: result ... structure containing the estimation results.
%         th     ... theta estimate
%         (Ac,Kc,Cc) ... state space system in canonical form.
%         Omegac ... noise variance estimate
%         thi   ... initial estimate corresponding to CCA.
%         thi2  ... initial estimate after RH-procedure.
%         lle   ... minimizing value of the unrestricted quasi scaled log
%                  likelihood (over M_n). 
%
% AUTHOR: dbauer, 27.10.2019

if nargin<8
    restrict.det_res = 0;
end
if nargin<7
    Pbull = -1;
end;

% calculate c as the total number of unit roots 
c = sum((urs(:,2)+1).*urs(:,3));

if c>n
    disp('Wrong specification: c cannot be larger than n!');
    result = est_result();
    return;
end;
% extract deterministics
dt = z(:,s+1:end);
m = size(dt,2);

y = z(:,1:s);

De = (dt\y)';
tily = y - dt*De';

if nargin<9
    % initial VAR estimate
    [k,sig,kbc,sigbc,phi]  = aicest(tily,size(y,2),nmax);
    
    % CVA estimate
    if 2*k*s<n
        k = n;
    end;
    
    [thi,Ai,Ki,Ci,Omegai] = CCA(tily,n,2*k,2*k,0);
    thi.D = De;
end;

% replace the estimate for the deterministics, if any 
if size(thi.D)~=m
    thi.D = De;
end


% find parameters in stationary version (in M_n)
parami = syst2param(thi,0,restrict);

% improve estimate
options = optimoptions('fminunc','display','final');
options.MaxFunctionEvaluations = 10000;
%restrict.scale = ones(length(parami),1);

pare = est_cal_like_hess(z,s,m,n,0,thi,restrict);

% random improvement of stable systems
result = compile_results(pare,s,m,n,0,restrict,z(:,1:s),z(:,s+1:end),Pbull,nmax);
lle_c = result.deviance;
lle=lle_c+1;
it = 0;
while ((it<5)&&(lle_c<lle))
    if it>0
        disp(sprintf('Running random improve %2d. Old value: %4.2f -> new value: %4.2f.',it,lle,lle_c));
    end
    lle = lle_c;
    result = random_improve_call(result,0.001,5);
    lle_c =result.deviance;
    it = it+1;
end;
pare = result.param;
lle =result.deviance; %cal_quasi_like(pare,z,s,m,n,0,Pbull,restrict);

[~,~,~,~,~,thi2] = param2syst(pare,s,n,m,0,restrict);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% now comes the MFI1 part.    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find next state space system with c common trends
[k]  = aicest(tily,size(y,2),nmax);

% find initial estimate.
[~,~,~,~,estimate] = VECM_MFI1(tily,k,S,urs,5);
[~,thi3] = RH_VECM_MFI1(z(:,1:s)-dt*thi2.D',S,thi2.A-thi2.K*thi2.C,thi2.K,estimate,5);

thi3.D = thi2.D;
thi3.B = zeros(n,m);
% find corresponding params
parami = th2param_MFI1(thi3,urs,1);
restrict.scale = ones(length(parami),1);

parc = est_cal_like_hess_MFI1(z,s,S,n,urs,thi3,restrict);
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
    %sizOm = s*(s+1)/2;
    %result.param = result.param((sizOm+1):end);
    llc_c =result.deviance;
    it = it+1;
end;

% write out results.
parc = result.param;
result = compile_results_MFI1(parc,s,S,n,urs,restrict,z(:,1:s),z(:,s+1:end),Pbull,nmax);
result.theta.ur='MFI(1)';
%[llc,resc] = cal_quasi_like(parc,z,s,m,n,c,Pbull,restrict);
%% provide corresponding estimate.
%[~,~,~,~,~,thc] = param2syst(parc,s,n,m,c,restrict);
thc = result.theta;
thc.B = zeros(n,m);
Ac = thc.A;
Kc = thc.K;
Cc = thc.C;
Omegac = thc.Omega;

% fill in everything into result structure
%result = est_result();
%result.c = c;
%result.n = n;
%result.s = s;
%result.y = y;
%result.dt = dt;
%result.theta = thc;
%result.deviance = llc;
%result.aic = llc + 2*length(param);
%result.bic = llc + log(size(resc,1))*length(param);
%result.res = resc;
%result.restrict = restrict;
%result.Pbull = Pbull;
%result.call = sprintf('SPECM(z,%d,%d,%d,%d,%d)',s,n,c,nmax,Pbull);
%result.param = param;



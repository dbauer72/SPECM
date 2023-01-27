function  [result,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM_I1(z,s,n,c,nmax,Pbull,restrict,thi);
% SPECM optimizes the Gaussian pseudo likelihood starting from an inital
% CCA estimated system.
% 
% SYNTAX: [result,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM(z,s,n,c,nmax,Pbull,restrict);
%
% INPUT:  z ... Tx(s+m) data: [y_t,d_t]
%         s ... integer; output dimension
%         n ... state order to use
%         c ... number of common trends
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

if nargin<7
    restrict.det_res = 0;
end
if nargin<6
    Pbull = -1;
end;

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
T = size(y,1); 

if nargin<8
    % initial VAR estimate
    [k,sig,kbc,sigbc,phi]  = aicest(tily,size(y,2),nmax);
    
    % CCA estimate
    if 2*k>sqrt(T) % reduce k if too big. 
        k = max(1,floor(sqrt(T)/2)); 
    end
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

% find parameters in stationary version
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

[~,~,~,~,~,the] = param2syst(pare,s,n,m,0,restrict);

% find next state space system with c common trends, if c>0. Else skip this
% step. 
if c>0 
    [LL,alphahat,betahat,Chat] = RH_specm(z(:,1:s)-dt*the.D',the.A-the.K*the.C,the.K,s-c,'y');

    Abar = the.A-the.K*the.C;
    thi2=the;
    the.A =Abar+the.K*Chat;
    the.C = Chat;

    % find initial estimate.
    % find corresponding params
    parami = syst2param(the,c,restrict);
    %restrict.scale = ones(length(parami),1);

    parc = est_cal_like_hess(z,s,m,n,c,the,restrict);
    param = parc;

    % random improvement of I(1) systems
else % if c==0 (stationary case)
    parc = pare;
end

% collect all results. 
result = compile_results(parc,s,m,n,c,restrict,z(:,1:s),z(:,s+1:end),Pbull,nmax);
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
result = compile_results(parc,s,m,n,c,restrict,z(:,1:s),z(:,s+1:end),Pbull,nmax);
result.theta.ur='I(1)';
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



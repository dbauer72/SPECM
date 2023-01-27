function  [result,thc,Ac,Kc,Cc,Omegac,thi,lle] = StSp_I0(z,s,n,nmax,Pbull,thi);
% StSp_I0 optimizes the Gaussian pseudo likelihood starting from an inital
% CCA estimated system.
% 
% SYNTAX: [result,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM(z,s,n,nmax,Pbull,thi);
%
% INPUT:  z ... Tx(s+m) data: [y_t,d_t]
%         s ... integer; output dimension
%         n ... state order to use
%         nmax ... maximum for AR order selection.
%         Pbull ... indicator; if Pbull>0, state is started in stationary
%                   distribution.
%         thi    ... use initial estimate rather than CCA.
%
% OUTPUT: result ... structure containing the estimation results.
%         th     ... theta estimate
%         (Ac,Kc,Cc) ... state space system in canonical form.
%         Omegac ... noise variance estimate
%         thi   ... initial estimate corresponding to CCA.
%         lle   ... minimizing value of the unrestricted quasi scaled log
%                  likelihood (over M_n). 
%
% AUTHOR: dbauer, 27.10.2019

restrict.det_res = 0;

if nargin<5
    Pbull = -1;
end;

% extract deterministics
dt = z(:,s+1:end);
m = size(dt,2);

y = z(:,1:s);

De = (dt\y)';
tily = y - dt*De';

if nargin<6
    % initial VAR estimate
    [k,sig,kbc,sigbc,phi]  = aicest(tily,size(y,2),nmax);
    
    % CCA estimate
    if 2*k*s<n
        k = n;
    end;
    
    [thi,Ai,Ki,Ci,Omegai] = CCA(tily,n,2*k,2*k,0);
    thi.D = De;
end;
% find parameters in stationary version
parami = syst2param(thi,0,restrict);

% improve estimate
options = optimoptions('fminunc','display','iter');
options.MaxFunctionEvaluations = 10000;
%restrict.scale = ones(length(parami),1);

%pare = est_cal_like_hess(z,s,m,n,0,thi,restrict);

% random improvement of stable systems
result = compile_results(parami,s,m,n,0,restrict,z(:,1:s),z(:,s+1:end),Pbull,nmax);
result = random_improve_call(result,0.001,5);
    
lle_c = result.deviance;
lle=lle_c+1;
it = 0;
while ((it<5)&&(lle_c<lle))
    if it>0
        disp(sprintf('Running random improve %2d. Old value: %4.2f -> new value: %4.2f.',it,lle,lle_c));
    end
    lle = lle_c;
    %result = random_improve_call(result,0.001,5);
    %lle_c =result.deviance;
    it = it+1;
end;
pare = result.param;
lle =result.deviance; %cal_quasi_like(pare,z,s,m,n,0,Pbull,restrict);

result.theta.ur='I(0)';
thc = result.theta;
thc.B = zeros(n,m);
Ac = thc.A;
Kc = thc.K;
Cc = thc.C;
Omegac = thc.Omega;


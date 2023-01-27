function  [result,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM_aom(z,s,n,c,nmax,Pbull,restrict);
% SPECM_aom optimizes the Gaussian pseudo likelihood starting from an inital
% CCA estimated system.
% 
% SYNTAX: [result,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM_aom(z,s,n,c,nmax,Pbull,restrict);
%
% INPUT:  z ... Tx(s+m) data: [y_t,d_t]
%         s ... integer; output dimension
%         n ... state order to use
%         c ... number of common trends
%         nmax ... maximum for AR order selection.
%         Pbull ... indicator; if Pbull>0, state is started in stationary
%                   distribution.
%         restrict ... structure specifying restrictions.
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
% REMARK: same as SPECM, but optimizes by alternating over th, Omega and D.
% 
% AUTHOR: dbauer, 27.10.2019


if nargin<7
    restrict.det_res = 0;
end
if nargin<6
    Pbull = 0;
end;

% extract deterministics
dt = z(:,s+1:end);
m = size(dt,2);

y = z(:,1:s);

De = (dt\y)';
tily = y - dt*De';

% initial VAR estimate
[k,sig,kbc,sigbc,phi]  = aicest(tily,size(y,2),nmax);

% CCA estimate
[thi,Ai,Ki,Ci,Omegai] = CCA(tily,n,2*k,2*k,1);
thi.D = De;
thi.Omega = Omegai; 
%
the = optim_quasi_like(thi,z,s,n,0,nmax,Pbull,restrict);
pare = syst2param(the,0,restrict);
[lle,rese] = cal_quasi_like(pare,z,s,m,n,0,Pbull,restrict);

% find parameters in stationary version
%parami = syst2param(thi,0,restrict);

% C does not contain parameters, A only in certain entries, B is full.
%Ae = thi.A;
%Ke = thi.K;
%Ce = thi.C;
%
%[n,s]=size(Ke);
%
%parami = th2param(thi,0);
%paraomi = extr_lowtri(thi.Omega);
%parami = [paraomi(:);parami(:);De(:)];

% improve estimate
options = optimoptions('fminunc','display','iter');
options.MaxFunctionEvaluations = 20000;
%restrict.scale = ones(length(parami),1);

%[pare,fval,exitflag] = fminunc(@(x) cal_quasi_like(x,z,s,m,n,0),parami,options);
%if exitflag>0
%    hess = hessian(@(x) cal_quasi_like(x,z,s,m,n,0),parc);
%    restrict.scale = 1./sqrt(diag(hess));
%    pare = pare.*sqrt(diag(hess));
%    [pare,fval,exitflag] = fminunc(@(x) cal_quasi_like(x,z,s,m,n,0),pare,options);
%end

%lle =cal_quasi_like(pare,z,s,m,n,0,0,restrict);
%
%[~,~,~,~,~,the] = param2syst(pare,s,n,m,0,restrict); 

%paro = pare(1:(s*(s+1)/2));
%Omegae = fill_lowtri(paro,s);
%pare(1:(s*(s+1)/2))=[];
%the = param2th(pare(1:(end-s*m)),s,n,0);
%the.Omega = Omegae;
%pare(1:(end-s*m)) = [];
%the.D = reshape(pare,s,m);

% find next state space system with c common trends
[LL,alphahat,betahat,Chat] = RH_vecm(z(:,1:s)-dt*the.D',the.A-the.K*the.C,the.K,s-c,'y');

Abar = the.A-the.K*the.C;
thi2=the;
the.A =Abar+the.K*Chat;
the.C = Chat;

% find initial estimate. 
% find corresponding params
thc = optim_quasi_like(the,z,s,n,c,nmax,Pbull,restrict);

%parami = syst2param(the,c,restrict);
%restrict.scale = ones(length(parami),1);

%paraomi = extr_lowtri(the.Omega);

% if restrict.det_res % restrict the highest order deterministic term.
%     the2 = param2th(parami,s,n,c);
%     C1 = the2.C(:,1:c);
%     D = the.D; 
%     pard = C1'*D(:,1);
%     D2 = D(:,2:end);
%     if m>1
%         pard = [pard(:);D2(:)];
%     end;
%     parami = [paraomi(:);parami(:);pard];
% else
%     parami = [paraomi(:);parami(:);the.D(:)];
% end
%th3 = param2th(parami(2:end),s,n,c)

% optimize full function
%[parc,~,exitflag] = fminunc(@(x) cal_quasi_like(x,z,s,m,n,c,Pbull,restrict),parami,options);
%if exitflag>0 
%    hess = hessian(@(x) cal_quasi_like(x,z,s,m,n,c,Pbull,restrict),parc);
%    restrict.scale = 1./sqrt(diag(hess));
%    parc = parc.*sqrt(diag(hess));
%    parc = fminunc(@(x) cal_quasi_like(x,z,s,m,n,c,Pbull,restrict),parc,options);
%end
%param = parc;

parc = syst2param(thc,c,restrict);
[llc,resc] = cal_quasi_like(parc,z,s,m,n,c,Pbull,restrict);
% provide corresponding estimate.
%[~,~,~,~,~,thc] = param2syst(parc,s,n,m,c,restrict); 
%paro = parc(1:(s*(s+1)/2));
%Omegac = fill_lowtri(paro,s);
%parc(1:(s*(s+1)/2))=[];
%
%nd = s*m;
%if restrict.det_res 
%    nd = nd- (s-c);
%end;
%
%thc = param2th(parc(1:(end-nd)),s,n,c);
%thc.Omega = Omegac;
%parc(1:(end-nd)) = [];
%
%if restrict.det_res
%    D = zeros(s,m);
%    pc = parc(1:c);
%    D(:,1)=thc.C(:,1:c)*pc(:);
%    parc(1:c)=[];
%    if m>1
%        D(:,2:end)= reshape(parc,s,m-1);
%    end
%else
%    D = reshape(parc,s,m);
%end
%thc.D = D; 
%%thc.LL = llc; 

Ac = thc.A;
Kc = thc.K;
Cc = thc.C;
Omegac = thc.Omega; 

% fill in everything into result structure
result = est_result();
result.c = c;
result.n = n;
result.s = s;
result.y = y;
result.dt = dt;
result.theta = thc;
result.deviance = llc;
result.res = resc;
result.restrict = restrict;
result.Pbull = Pbull;
result.call = sprintf('SPECM_aom(z,%d,%d,%d,%d,%d)',s,n,c,nmax,Pbull);
result.param = parc;



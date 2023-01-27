function  [result,thc,Ac,Kc,Cc,Omegac] = SPECM_I2(z,s,n,ds,nmax,Pbull,restrict,thi);
% SPECM optimizes the Gaussian pseudo likelihood starting from an inital
% CCA estimated system.
% 
% SYNTAX: [result,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM_I2(z,s,n,ds,nmax,Pbull,restrict);
%
% INPUT:  z ... Tx(s+m) data: [y_t,d_t]
%         s ... integer; output dimension
%         n ... state order to use
%         ds ... 2x1 number of [I(2),I(1)] common trends
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
% Correspondence indices of Paruolo and Johansen rs and state space unit
% root structure ds:
%    rs = [rs1,rs2] where rs1 is the rank of Phi, rs2 the rank of
%    alpha_perp Gamma beta_perp. Thus rs1 is the number of stationary
%    components in y_{t-1}, if Delta y_{t-1} is regressed out, s- rs1 - rs2 the
%    number of I(2) common trends. 
%  ds = [ds1,ds2] where ds1 is the number of cointegrating relations, ds2
%  the number of additional I(1) trends not contained in Delta y_t. 
%
% Thus: ds1 = s-rs1-rs2, ds2 = rs2. 
%  Or: rs1 = s-ds1 -ds2, rs2 = ds2. 
%
% AUTHOR: dbauer, 27.10.2019

if length(ds)==1 % spec might be for I(1) -> add one I(2) common trend. 
    ds = [1,ds];
end;


if sum(ds)>s
    error(' SPECM_I2: Too many common trends, system would not be minimal!');
end;

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

%%%%%%%%%%%%%%%%%%%%%
% initial  estimate %
%%%%%%%%%%%%%%%%%%%%%
% short variant: 
%thi = initial_I2(y,rs,n,nmax);
if nargin<8 
    % long one: 
    [k,sig,kbc,sigbc,phi]  = aicest(tily,size(y,2),nmax);
    %
    %% CCA estimate
    if k<n
        k = 2*n;
    end;

    rs(1)=s-sum(ds);
    rs(2)=ds(2);
    c = 2*ds(1)+ds(2); 

    if n< c
        n = c;
        disp(' SPECM_I2: too small system order -> correcting.');
    end;


    [thi,Ai,Ki,Ci,Omegai] = CCA(tily,n,2*k,2*k,0);
end

% replace the estimate for the deterministics, if any 
if size(thi.D)~=m
    thi.D = De;
end

% find corresponding params for stationary system
parami = syst2param(thi,0,restrict);
restrict.scale = ones(length(parami),1);
[lli,resi] = cal_quasi_like(parami,y,s,m,n,0,Pbull,restrict);
thi.Omega = resi'*resi/size(y,1); 
parami(1:s*(s+1)/2) = extr_lowtri(thi.Omega);

restrict.scale = ones(length(parami),1);
options = optimoptions('fminunc','display','iter','MaxIterations',1000);
options.MaxFunctionEvaluations = 20000;
[parami,~,~] = fminunc(@(x) cal_quasi_like(x,tily,s,m,n,0,0,restrict),parami,options);

[lli,resi] = cal_quasi_like(parami,tily,s,m,n,0,Pbull,restrict);
thi = param2th(parami(s*(s+1)/2+1:end),s,n,0,restrict);
thi.Omega = resi'*resi/size(y,1);

%% approximate by I(2) system: 
%Abar = thi.A - thi.K*thi.C;
%B = thi.K;
%thi =  RH_VECM_LM(y,Abar,B,rs,0);
%thi.D = De;

thi = initial_I2(tily,rs,n,nmax);
thi.D = De;
thi.m = 1;
thi.B = zeros(n,1);

% find paramters 
parami = syst2param_I2(thi,ds,1);
restrict.scale = ones(length(parami),1);

% improve estimate
%options = optimoptions('fminunc','display','iter');
options = optimoptions('fminunc','display','iter','MaxIterations',1000);
options.MaxFunctionEvaluations = 20000;

[parc,fval,exitflag] = fminunc(@(x) cal_quasi_like_I2(x,y,s,m,n,ds,0,restrict),parami,options);

[llc,resc] = cal_quasi_like_I2(parc,y,s,m,n,ds,Pbull,restrict);
% provide corresponding estimate.
[~,~,~,~,~,thc] = param2syst_I2(parc,s,n,m,ds,restrict); 
thc.B = zeros(n,m);

Ac = thc.A;
Kc = thc.K;
Cc = thc.C;
Omegac = thc.Omega; 

thc.ur = 'I(2)';
thc.urs = ds;

% fill in everything into result structure
result = est_result();
result.urs = ds;
result.ur = 'I(2)';
result.n = n;
result.s = s;
result.y = y;
result.dt = dt;
result.theta = thc;
result.deviance = llc;
result.aic = llc + 2*length(parc);
result.bic = llc + log(size(resc,1))*length(parc);
result.res = resc;
result.restrict = restrict;
result.Pbull = Pbull;
result.call = sprintf('SPECM_I2(z,%d,%d,%d,%d,%d)',s,n,ds(2),nmax,Pbull);
result.param = parc;



function the = optim_quasi_like(thi,z,s,n,c,nmax,Pbull,restrict);
% optim_qasi_like implements the alternating optimization algorithm where
% iteratively the quasi likelihood is optimized over Omega, D and th. 
%
% SYNTAX: the = optim_quasi_like(thi,z,s,n,c,nmax,Pbull,restrict);
%
% INPUT: thi ... initial estimate in theta format
%        z   ... Tx(s+m) observations.
%        s   ... integer; dimension of endogenous variable.
%        n   ... integer; state dimension
%        c   ... integer; number of common trends.
%        nmax ... integer; maximal number of states (only needed for
%                consistency of input params)
%        Pbull ... indicator; if Pbull>0 then state is started from
%                  stationary distribution.
%        restrict ... structure spefifying restrictions.
%
% OUTPUT: the ... theta structure of estimated system
%
% AUTHOR: dbauer, 14.1.2020.

dt = z(:,s+1:end);
m = size(dt,2);


options = optimoptions('fminunc','display','iter');
options.MaxFunctionEvaluations = 20000;
%restrict.scale = ones(length(parami),1);

the = thi; 
Ome = thi.Omega; 
paraom = extr_lowtri(thi.Omega); % parameters for Omega.
sizOm = length(paraom); 

De = thi.D;
parath = th2param(thi,c,0); % parameters for theta structure. 

% extract parameters for D 
if restrict.det_res == 0
    nparD = s*m;
else
    nparD = s*(m-1)+c;
end
parde = thi.D(:); %parath(end-nparD+1:end);
%parath(end-nparD+1:end)=[];

% iterate cycles over three elements for 5 iterations (hard coded!!) 
for k=1:5
    % extract parameters.
    parath = th2param(the,c,1);
    % optimize over th
    [pare,fval,exitflag] = fminunc(@(x)  cal_quasi_like_th(x,the,z,s,m,n,c,Pbull,restrict),parath,options);
    the = param2th([pare(:);parde(:)],s,n,c,restrict);
    the.Omega = Ome;
    
    % optimize over D
    if (restrict.det_res)&&(c>0) % restrict the highest order deterministic term.    
        C1 = the.C(:,1:c);
        D = De; 
        parde = C1'*D(:,1);
        D2 = D(:,2:end);
        if size(D,2)>1
            parde = [parde(:);D2(:)];
        end;
    else
        parde = De(:);
    end
    [parde,fval,exitflag] =  fminunc(@(x) cal_quasi_like_D(x,the,z,s,m,n,c,Pbull,restrict),parde,options);
    De =  fill_D(parde,the.C(:,1:c),s,m,restrict);
    the.D = De;
    % optimize over Omega.
    paraom = extr_lowtri(the.Omega);
    [pare,fval,exitflag] = fminunc(@(x)  cal_quasi_like_Omega(x,the,z,s,m,n,c,Pbull,restrict),paraom,options);
    Ome =  fill_lowtri(pare,s);
    the.Omega = Ome;
end; 

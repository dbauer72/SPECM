function result = random_improve_call(result,k,M);
% random_improve_call perturbs the parameter vector in result 
% using a normally distributed disturbance with variance k
% and optimizes the criterion function.
% M random draws are chosen. 
%
% INPUT: result ... result structure
%
% OUTPUT: result ... updated result structure
%
% AUTHOR: dbauer, 27.3.2020.

z = [result.y,result.dt];
s=result.s;
m = size(result.dt,2);
n=result.n;
urs = result.urs;
Pbull=result.Pbull;
rest = result.restrict;


% evaluate initial
param = result.param;
if Pbull<0
    sizOm = s*(s+1)/2;
    param = param(sizOm+1:end);
end;
param_old = param;

ur = result.ur; 
switch ur
    case 'I(1)' 
        c = urs;
        ll_cur = cal_quasi_like(param,z,s,m,n,c,Pbull,rest);
    case 'I(2)' 
        c = urs; 
        ll_cur = cal_quasi_like_I2(param,z,s,m,n,c,Pbull,rest);
    case 'MFI(1)'
        param = result.param;
        param_old = param;
        urs = result.urs;
        S = result.S;
        ll_cur = cal_quasi_like_MFI1(param,z,s,S,n,urs,Pbull,rest);
    otherwise
        ll_cur = cal_quasi_like(param,z,s,m,n,0,Pbull,rest);
end

options = optimoptions('fminunc','display','iter','MaxIterations',1000);
options.MaxFunctionEvaluations = 20000;
    
h = waitbar(0,'Please wait ...','Name','Random_improve');

for j=1:M
    waitbar(j/M,h,sprintf('%d of %d',j,M));
    try
        ll_n= ll_cur;
    switch ur 
        case 'I(1)' 
            [parn,exitflag] = fminunc(@(x) cal_quasi_like(x,z,s,m,n,c,Pbull,rest),param+randn(length(param),1)*k,options);
            ll_n = cal_quasi_like(parn,z,s,m,n,c,Pbull,rest);
        case 'I(2)'
            [parn,exitflag] = fminunc(@(x) cal_quasi_like_I2(x,z,s,m,n,c,Pbull,rest),param+randn(length(param),1)*k,options);
            ll_n = cal_quasi_like_I2(parn,z,s,m,n,c,Pbull,rest)
        case 'MFI(1)'
            [parn,exitflag] = fminunc(@(x) cal_quasi_like_MFI1(x,z,s,S,n,urs,Pbull,rest),param+randn(length(param),1)*k,options);
            ll_n = cal_quasi_like_MFI1(parn,z,s,S,n,urs,Pbull,rest);
        otherwise
            [parn,exitflag] = fminunc(@(x) cal_quasi_like(x,z,s,m,n,0,Pbull,rest),param+randn(length(param),1)*k,options);
            ll_n = cal_quasi_like(parn,z,s,m,n,0,Pbull,rest);
    end
    end
    if ll_n< ll_cur
        param = parn;
        ll_cur = ll_n;
    end
end
close(h); 

% update result object. 
if norm(param-param_old)>10^(-10) % changes have been achieved
    if Pbull<0
        n_par_om = 0;
    else 
        n_par_om = s*(s+1)/2;
    end     
    switch ur
        case 'I(1)'
            [th,par] = param2th(param(n_par_om+1:end),s,n,c,rest);
            [llc,resc] =  cal_quasi_like(param,z,s,m,n,c,Pbull,rest);
        case 'I(2)'
            [th,par]  = param2th_I2(param(n_par_om+1:end),s,n,c,rest);
            [llc,resc] =  cal_quasi_like_I2(param,z,s,m,n,c,Pbull,rest);
        case 'MFI(1)'
            [th,par]  = param2th_MFI1(param(1:end),s,n,urs,rest);
            [llc,resc] =  cal_quasi_like_MFI1(param,z,s,S,n,urs,Pbull,rest);
        otherwise
            [th,par]  = param2th(param(n_par_om+1:end),s,n,0,rest);
            [llc,resc] =  cal_quasi_like(param,z,s,m,n,0,Pbull,rest);
    end
    th.Omega = resc'*resc/size(z,1);

    th.B = zeros(n,m);
    if isempty(par)
        th.D = zeros(s,m);
    else
        if (length(par)==s*m)
            th.D = reshape(par,s,m);
        else
            th.D = result.theta.D;
        end
    end
    result.theta= th;
    result.deviance = llc;
    result.res = resc;
    switch ur 
        case {'I(1)','I(2)','I(0)'}
            if (Pbull<0)
                param = [extr_lowtri(th.Omega);param];
         end;
    end
    result.param = param;
    result.aic = llc+2*length(param);
    result.bic = llc+log(size(z,1))*length(param);
end



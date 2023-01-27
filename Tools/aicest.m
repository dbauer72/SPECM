function [k,sig,kbc,sigbc,phi,sigphi,thar] = aicest(z,s,nmax);
%
% 	estimates the (column) truncation index using AIC
%
% SYNTAX: [k,sig,kbc,sigbc,phi,sigphi,thar] = aicest(z,s,nmax);
%
% INPUTS:   z ... Tx (s+m); [y u], observations output (y) and input (u)
%		    s ... dim(y)
%		  nmax .. maximal lag length
%
% OUTPUTS:  k   ...  estimated lag order according to AIC
%		    sig ...  AIC values
%           kbc ... lag order estimate using BIC
%         sigbc ... criterion values for BIC
%           phi ... AR estimate using lag length k. 
%        sigphi ... scaled log likelihood values for all lag orders. 
%          thar ... theta structure corresponding to AR(k).
%		  
%  AUTHOR: adaptations by dbauer, 13.1.2020.

[T,nz] = size(z);
nu = nz - s;

if (nmax<0)
    k_est = -nmax;
    nmax = -nmax;
else
    k_est = [];
end

Teff = T-nmax; % effective sample size.
Xfull = zeros(Teff,nmax*nz);
for i=1:nmax
    Xfull(:,(i-1)*nz+[1:nz]) = z(nmax-i+[1:Teff],:);
end
y = z(nmax+[1:Teff],1:s);

sig(1) = log(det(y'*y/Teff));
sigphi(1) = sig(1);
sigbc(1) = sig(1);

for i=1:nmax
    res_i = y - Xfull(:,1:nz*i)*(Xfull(:,1:nz*i)\y);
    Om_i = res_i'*res_i/Teff;
    dOm_i = det(Om_i);
    if dOm_i<10^(-20)
        dOm_i = 10^(-20);
    end;
	sig(i+1) = log(dOm_i) + 2*i*s*nz/Teff;
	sigphi(i+1) = log(dOm_i);
	sigbc(i+1) = log(dOm_i) + i*s*nz*log(Teff)/Teff;
end;

% AIC/BIC defined as T*log det hat Sigma + T*s + C_T num_par  
% rescale 
sig = sig*Teff+Teff*s;
sigbc = sigbc*Teff+Teff*s; 

k = find(~(sig > min(sig)))-1;
kbc = find(~(sigbc > min(sigbc)))-1;
if isempty(k_est)
    k_est = k;
end;

phi = (Xfull(:,1:nz*k_est)\y)';
if k_est>0
    if nargout>6 
    thar = theta_urs();
    thar.which = 'poly'; %use polynomial form.
    % find indices for y and z terms
    indy = [1:s];
    if nu>0
        indz = [(s+1):nz];
    else
        indz = [];
    end
    for j=2:k_est
        indy = [1:s,indy+nz]; 
        if nu>0
            indz = [(s+1):nz,indz+nz];
        end
    end;
    thar.a = [eye(s),-phi(:,indy)];
    thar.d = phi(:,indz);
    thar.b = eye(s);
    res_i = y - Xfull(:,1:nz*k_est)*phi';
    thar.Omega = res_i'*res_i/Teff;
    thar.m = nu;
    thar.num_param = length(phi(:));
end
else
    thar.a = eye(s);
    thar.d = zeros(s,nu);
    thar.b = eye(s);
    thar.Omega = y'*y/Teff;
    thar.m = nu;
    thar.num_param = 0;
end



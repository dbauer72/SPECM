function [th,A,K,C,Omega] = CCA(z,n,kcol,krow,plots);
% CCA implements canonical correlations analysis subspace method.
%
% SYNTAX: [th,a,b,c,d,k,Omega] = CCA(z,n,kcol,krow,p);
%
% INPUTS:   z ... Tx s matrix of observations.
%           n ... integer; system order; if n=[] estimated using SVC.
%           kcol, krow ... integers; past and future horizons for CCA.
%           plots ... indicator; if plots>0 singular values are plotted. 
%
% OUTPUTS:  th ... theta structure of estimated system.
%           (A,K,C) ... estimated state space system
%           Omega ... sxs estimated innovation variance matrix.
% 
% REMARK: CCA subspace algorithm (see Larimore 1987), no exogenous inputs; uses a regression framework
%
% dbauer, 27.10.2019


[T,nz]=size(z);

% ------ data Hankel matrices -----
while krow*nz> (T-krow-kcol)
    krow = krow-1;
end

for i=1:krow
   Yf((i-1)*nz+[1:nz],:) = z(kcol+i+[0:T-krow-kcol],1:nz)';
end;
for i=1:kcol
   Zp((kcol-i)*nz+[1:nz],:) = z(i+[0:T-krow-kcol],1:nz)';
end;
      
% ------ regressions to obtain hat beta ----
betaz = Yf/Zp;

% ---- SVD ----
Wf2 = (Yf*Yf')/T;
Wf = inv(chol(Wf2)');

Wp2 = (Zp*Zp')/T;
if min(eig(Wp2))<10^(-6)
    Wp2 = Wp2 + eye(size(Wp2))*10^(-6);
end
Wp =  chol(Wp2)';

[U,S,V] = svd(Wf*betaz*Wp);

% order estimation using SVC
if isempty(n)|(plots)
    s = diag(S);
    svc = s.^2 + log(T)/T*2*[1:length(s)]'*nz;
    nmax = length(s)-1;
    [minn,nhat] = min(svc(1:nmax+1));
    if (nhat>1)
       nhat = nhat-1;
    else
        nhat=1;
    end;
    
    nhat = max(nhat,1); % make sure that the estimated order is not too small.
end

if plots
    
    figure;
    plot(s,'x');
        hold on;
    plot([0,nmax],sqrt(1-log(T)/T)*[1,1],'m');
    chat = sum(s>sqrt(1-log(T)/T));
    title(sprintf('Singular values: c: %d, SVC: %d',chat,nhat));
    plot(nhat,S(nhat,nhat),'ro');
    for c=1:chat,
        plot(c,S(c,c),'mo');
    end;
end;

if isempty(n)
    n = nhat;
end;
beta = V(:,1:n)'*inv(Wp);
beta = beta(:,1:n)\beta;

% ------ system matrices ------
x = beta*Zp;
C = Yf(1:nz,:)/x;
res = Yf(1:nz,:) - C*x;
AK = x(:,2:end)/[x(:,1:end-1);res(:,1:end-1)];
A = AK(:,1:n);
K = AK(:,n+1:end);
Omega = (res*res')/T;

% ----------   transformation to Echelon form
% generic neighbourhood
th = ss2ech_n(A,K,C);
th.Omega = Omega; 

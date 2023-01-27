function [estimate,u,Dhat] = update_estimate(estimate,pX,pXk,Xm2,Joh_j);
% final regression in Johansen Schaumburg proc fixing betas and estimating
% jointly alphas as well as Gammas.
%
% SYNTAX: estimate = update_estimate(estimate,pX,pXk,Xm2,m);
%
% INPUT: estimate ... vecm_mfi object. 
%        pX   ... T x s dependent var p(L)X_t
%        pXk  ... T x sk regressors p(L)X_{t-j}
%        Xm2  ... T x s x M regressors X_t^(m)
%     
% OUTPUT: estimate ... updates estimate object.
%         u        ... T x n_u real matrix of deterministic regressors.
%         Dhat     ... s x n_u real matrix; corresponding estimates. 
%
% AUTHOR: dbauer, 26.11.2015. 

[T,s]= size(pX);
M = size(Xm2,3);
ind = estimate.k*s;
if Joh_j<4
    ind = ind+1;
    u = pXk(:,ind);
else
    u = zeros(T,0);
end;
X = pXk(:,1:ind);

% generate regressors
for m=1:M
    xx = squeeze(Xm2(:,:,m));
    if (m>1)&&(m<M)  % complex root 
        X = [X,xx*estimate.beta{m}];        
    else 
        beta = estimate.beta{m};
        if size(beta,1)>s 
            xxh = xx(:,[1:s,2*s+1])*estimate.beta{m};
        else
            xxh = xx(:,1:s)*estimate.beta{m};
        end
        X = [X,xxh];       
    end;
    ind = [ind,size(X,2)];    
end;

% OLS
bhat = (X\pX)';

% fill in results in object.
res = pX - X * bhat';
estimate.Omega = res'*res/T;
if Joh_j<4
    estimate.Gamma = bhat(:,1:(ind(1)-1));
else
    estimate.Gamma = bhat(:,1:ind(1));
end

for m=1:M
    estimate.alpha{m} = bhat(:,[(ind(m)+1):ind(m+1)]);
end;

% now deal with deterministics
Dhat = zeros(s,0);
if Joh_j<4
    Dhat = bhat(:,ind(1));
end;

if Joh_j==4
    Dhat(:,end+1) = estimate.alpha{1}*estimate.beta{1}(s+1,:)';
    u = ones(T,1);
end

if Joh_j<5
    % complex roots.
    for m=2:(M-1)
        xx = squeeze(Xm2(:,:,m));
        u = [u,xx(:,[s+1,end]);];
        alpha = estimate.alpha{m};
        beta = estimate.beta{m};
        
        Dhat = [Dhat,alpha*beta([s+1,end],:)'];
    end
    
    % root at z=-1.
    xx = squeeze(Xm2(:,:,M));
    u = [u,xx(:,2*s+1)];
    Dhat = [Dhat,estimate.alpha{M}*estimate.beta{M}(s+1,:)'];
end

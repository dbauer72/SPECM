function [pX,pXk,Xm,fr,urs] = cal_reg_urs(y,S,k); 
% provides the regressors X_t^{(m)} as well as p(L)X_t 
% according to Johansen and Schaumburg.

% SYNTAX: [pXk,Xm,fr] = cal_reg_urs(y,urs,M,k,S);
%
% INPUT:  y   ... Txs real time series.
%         S   ... integer; number of seasons. 
%         k   ... max lag order.
%         Joh_j ... specification according to Johansen (1997).
%
% OUTPUT: pX ... T x s: p(L)X_t
%         Xm  ... T x s x S array of regressors X_t^{(m)}
%         fr  ... vector of frequencies. 
%         urs ... initial unit root structure without number of common
%                cycles.
%
% REMARK: determinstics: Joh_j: 1,2,3 ... unrestricted constant, seasonal
% terms added to cointegrating relations, where unit roots are present,
% omitted else. 
%   Joh_j: 4 ... also for z=1 if unit root present cointegrating vector
%   adapted. 
%  Joh_j: 5 ... no determinstics included. 
%
% AUTHOR: dbauer, 26.11.2015, adapted from original from Martin Wagner.
% updated for general S: 19.8.2020.


% calc frequencies.
[T,s] = size(y);
om = [0:(S/2)]'/S*2; % frequencies between [0,1]*pi.
fr = exp(pi*om*sqrt(-1)); % corr. unit roots. 

urs(:,1)= om(:);
urs(:,2)= 1;
urs([1,end],2)=0;

% filter regressors
pX = y;
for f=1:length(fr)
    Xm2(:,:,f)=y(1:end-1,:)/fr(f); % X(t-1)/z_m.
end;

for f=1:length(fr)
    pX = filter([1,-conj(fr(f))],1,pX);
    if urs(f,2) % for complex frequency -> also filter conjugate root.
        pX = filter([1,-(fr(f))],1,pX);
    end;
    for g=1:length(fr)
        if (g~=f)
            Xm2(:,:,g) = filter([1,-conj(fr(f))],1,Xm2(:,:,g))/(1-conj(fr(f))*fr(g));
            if urs(f,2)
                Xm2(:,:,g) = filter([1,-(fr(f))],1,Xm2(:,:,g))/(1-fr(f)*fr(g));
            end;
        else 
            if urs(g,2)
                Xm2(:,:,g) = filter([1,-(fr(f))],1,Xm2(:,:,g))/(1-(fr(f)*fr(g)));
            end;
        end;
    end;    
end;

pX = real(pX);

% transform to real ones and lag once. 
for j=1:size(Xm2,3)
    Xcomp = squeeze(Xm2(:,:,j));
    if urs(j,2)== 1 % imag root 
        Xm(:,1:2*s,j)=[NaN*ones(1,2*s);[real(Xcomp),imag(Xcomp)]];
    else 
        Xm(:,1:2*s,j)=[NaN*ones(1,2*s);[real(Xcomp),imag(Xcomp)]];
    end;
end;

% lag p(L)X up to maximum of k lags.
pXk = mlag(pX,k,NaN);
ii = 1:(s*k);
II = reshape(ii,k,s)';
pXk = pXk(:,II(:));

% eliminate rows with NaNs.
in = find(sum(isnan([pX,pXk]),2)==0);
in = setdiff(in,1:S); % take out first year if not done automatically. 

pX = pX(in,:);
pXk = pXk(in,:);
for j=1:size(Xm,3)
    XX = squeeze(Xm(:,:,j));
    XX = XX(in,:);
    Xm3(:,:,j)=XX;
end;

Xm = Xm3;



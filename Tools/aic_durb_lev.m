function [alpha,sigma,res,aic] = aic_durb_lev(y,lag,plots);
% pacf_sample calculates the pacf and plots it.
%
% SYNTAX: pacf_sample(y,lag,plots);
%
% INPUT:  y ... Txs real observation matrix.
%         lag ... integer; maximal number of lag
%
% OUTPUT: alphaest ... lagx1 estimated PACF.
%
% REMARK: uses the Durbin-Levinson iterative algorithm.
%
% AUTHOR: dbauer, 3.4.2020.

[T,nz] = size(y);
if (nz>1)
    disp(' Only SISO case!');
    return;
end;

R = mcovf(y,lag+1);
phip = 1;
alphaest(1)=1;
sig(1) = R(1);
for i=2:lag+1
    alphaest(i) = -phip*R(i:-1:2)'/sig(i-1);
    sig(i)=sig(i-1)*(1-alphaest(i)^2);
    phip = [1,phip(2:end)+phip(end:-1:2)*alphaest(i),alphaest(i)];
end;

alpha = phip;
res = y;
for j=1:lag
    res(lag+1:end) = res(lag+1:end) - alpha(j)*y((lag+1-j):end-j);
end
res(1:lag)= NaN;
sigma = res(lag+1:end)'*res(lag+1:end)/(T-lag);
aic = log(sig(:))+ 2*[0:lag]'/T;
if plots
    figure;
    plot([0:lag],sig,'x');
    hold on;
    plot([0:lag],aic,'r');
    [~,mi] = min(aic);
    plot(mi-1,aic(mi),'k*');
end

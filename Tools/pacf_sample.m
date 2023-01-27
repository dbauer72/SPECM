function alphaest = pacf_sample(y,lag,lab,plots,yhat,time);
% pacf_sample calculates the acf and the pacf and plots it.
%
% SYNTAX: pacf_sample(y,lag,plots);
%
% INPUT:  y ... Txs real observation matrix.
%         lag ... integer; maximal number of lag
%         lab ... string; name of series
%         plots ... indicate whether a plot is wanted.
%         yhat ... Tx1 real; extracted trend terms. 
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
end

if nargin< 6
    time = 1:T;
end

R = mcovf(y,lag+1);
phip = 1;
sig(1) = R(1);
for i=2:lag+1
    alphaest(i-1) = -phip*R(i:-1:2)'/sig(i-1);
    sig(i)=sig(i-1)*(1-alphaest(i-1)^2);
    phip = [1,phip(2:end)+phip(end:-1:2)*alphaest(i-1),alphaest(i-1)];
end;

acf = R(1:end)/R(1);

if plots 
    % plot y and extracted trend terms
    subplot(2,2,1);
    plot(time,y+yhat);
    hold on;
    plot(time,yhat);
    title(lab);
    hold off;
    
    % plot only y
    subplot(2,2,2);
    plot(time,y);
    title(lab);

    % plot acf
    subplot(2,2,3);
    
    % acf 
    stem([0:lag],acf);
    title(sprintf('ACF: %s',lab));

    hold on;
    plot([-1,lag+1],2/sqrt(T)*[1,1],'b--');
    plot([-1,lag+1],-2/sqrt(T)*[1,1],'b--');
    set(gca,'xlim',[-1,lag+1],'ylim',[-1.1,1.1]);    
    hold off;

    % pacf
    subplot(2,2,4);
    stem([1:lag],-alphaest);
    title(sprintf('PACF: %s',lab));

    hold on;
    plot([0,lag+1],2/sqrt(T)*[1,1],'b--');
    plot([0,lag+1],-2/sqrt(T)*[1,1],'b--');
    set(gca,'xlim',[0,lag+1],'ylim',[-1.1,1.1]);    
    hold off;
end
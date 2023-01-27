function ok = pres_select(which);
% presents the results from the selection process (callback to analyze_GUI).
%
% INPUT:  which ... either 'stat' or 'non-stat'. 
%
% AUTHOR: dbauerm 28.8.2020.

g = get(gcf,'userdata');
y = g{1};
dt = g{2};
time = g{4};


if strcmp(which,'non-stat')
    gg= findobj('tag','results');
    g = get(gg,'userdata');
    ord_obj = findobj('tag','Order_stat');
    ord = str2num(get(ord_obj,'string'));
else 
    gg= findobj('tag','results_stat');
    g = get(gg,'userdata');
    ord_obj = findobj('tag','Order');
    ord = str2num(get(ord_obj,'string'));
end

if ~isempty(g)
    results = g; 
    for j=1:length(results)
        sig(j)=results(j).aic;
        sigbc(j) = results(j).bic;
    end;
    
    figure(3)
    clf;

    plot([0:length(sig)-1],sig,'r');
    hold on;
    plot([0:length(sig)-1],sig,'rx');
    [~,k] = min(sig);
    plot([0:length(sig)-1],sigbc,'b');
    plot([0:length(sig)-1],sigbc,'bx');
    plot(k-1,sig(k),'ro');
    [~,kb] = min(sigbc);
    plot(kb-1,sigbc(kb),'bo');

    legend('AIC','BIC','bAIC','bBIC');
    title('AIC/BIC over order of state space system (log cond. likelihood)');

else
    msgbox('No results contained currently!','pres_select');
end;

try
    figure(1);    
    th_ss = results(ord+1).theta;
    dt = results(ord+1).dt;
    present_est_results(y,dt,time,th_ss);
end
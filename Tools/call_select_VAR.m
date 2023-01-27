function thar = call_select_VAR;
% script for selecting the lag length in the VAR case for analyze_GUI.
%
% OUTPUT:   thar ... theta structure of estimate.
%
% AUTHOR: dbauer, 18.8.2020.

g = get(gcf,'userdata');
y = g{1};
dt = g{2};
S = g{3};
time = g{4};
labs = g{5}; 

try
    e0 = findobj(gcf,'tag','T0');
    e1 = findobj(gcf,'tag','T1');
    ve0 = get(e0,'value');
    ve1 = get(e1,'value');
    if S>1
        es = findobj(gcf,'tag','TS');
        ves = get(es,'value');
    end    
catch
    disp('Wrong entry in field detrend.');
    return;
end

% detrend before estimating the model. 
T = size(y,1);
fdt = dt;
if ve0>0 
    fdt = [fdt,ones(T,1)];
end
if ve1>0 
    fdt = [fdt,[1:T]'];
end

if ((ves>0)&&(S>1))
    if (ve0>0)
        st = 2;
    else
        st = 1;
    end;
    for j=st:S
        dums = zeros(T,1);
        dums([j:S:T])=1;
        fdt = [fdt,dums];
    end
end

if ~isempty(dt)
    fdt = dt(end-length(y)+1:end,:);
else
    fdt = zeros(length(y),0);
end

% protect from double usage of deterministics. 
[q,r]= qr(fdt);
dr = diag(r);
cc = sum(dr>10^(-12));
fdt = q(:,1:cc);

yhat = fdt*(fdt\y);
ycur = y- yhat;


% now find maximum lag 
gg = findobj(gcf,'tag','MaxL');
nmax = str2num(get(gg,'String'));

% estimate using AIC 
[k,sig,kbc,sigbc,phi,sigphi,thar] = aicest(ycur,size(ycur,2),nmax);

% put in result. 
gg = findobj(gcf,'tag','Lag');
set(gg,'String',num2str(k));

figure(3)
clf;

plot([0:length(sig)-1],sig,'r');
hold on;
plot([0:length(sig)-1],sigbc,'b');
plot(k,sig(k+1),'ro');
[~,kb] = min(sigbc);
plot(kb-1,sigbc(kb),'bo');

legend('AIC','BIC','bAIC','bBIC');
title('AIC/BIC over VAR lag length (YW estimates)');

present_est_results(ycur,zeros(T,0),time,thar,0,labs);



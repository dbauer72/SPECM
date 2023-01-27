function th_ss = call_select_StSp;
% script for selecting the lag length in the VAR case for analyze_GUI.
%
% OUTPUT:   th_ss ... theta structure of estimate.
%
% AUTHOR: dbauer, 18.1.2023.

g = get(gcf,'userdata');
y = g{1};
dt = g{2};
S = g{3};
time = g{4};

sdet = 0; 
try
    odet = findobj(gcf,'tag','det_stat');
    sdet = get(odet,'value');
catch
    disp('Wrong entry in field deterministic.');
    return;
end


% get specification for selection
gg = findobj(gcf,'tag','MaxO_stat');
nmax = str2num(get(gg,'String'));

gg = findobj(gcf,'tag','Order_stat');
n = str2num(get(gg,'String'));

gg = findobj(gcf,'tag','crit_stat');
cr =get(gg,'Value');

gg = findobj(gcf,'tag','PEM_stat');
Pbull = 1-get(gg,'Value');


% include deterministics 
[T,s]= size(y);


% add deterministics to dt
if sdet
    t0o = findobj(gcf,'tag','T0');
    t0 = get(t0o,'value');
    if t0>0
        dt = [dt,ones(T,1)];
    end
    t1o = findobj(gcf,'tag','T1');
    t1 = get(t1o,'value');
    if t1>0
        dt = [dt,[1:T]'/T];
    end
    t2o = findobj(gcf,'tag','TS');
    t2 = get(t2o,'value');
    if t2>0
        ti = [1:T]';
        for j=0:S-1
            dt = [dt,(rem(ti,S)==j)];
        end
    end
end


results(1) = null_model(y,dt);
results(1).time =time;

h=  waitbar(0,'Iterating over order ...','Name','Order estimation');
for n=1:nmax
    waitbar(n/nmax,h,sprintf('Estimating order %d of %d.',n,nmax));
    [results(n+1)] = StSp_I0([y,dt],s,n,floor(sqrt(T)),Pbull);
    results(n+1).time =time;
    result = results(n+1);
    if results(n).deviance<result.deviance % system of order one less is better!
        result = results(n);
        th_est = result.theta;
                
        th_est.A = [th_est.A,zeros(n-1,1);zeros(1,n-1),0.1];
        th_est.K(end+1,:) = randn(1,s);
        th_est.C(:,end+1)=randn(s,1)*0.0001;
        th_est.B(end+1,:)= zeros(1,size(dt,2));
        result.theta= th_est;
        result.n = n;
        par = th2param(result.theta,0,1);
        s2 = s*(s+1)/2;
        result.param = [results(n).param(1:s2);par;result.theta.D(:)];
        results(n+1) = random_improve_call(result,.01,5);
     end
end
close(h);
        
gg= findobj('tag','results_stat');
set(gg,'userdata',results);

pres_select('stat');

% evaluate best model 
for j=1:length(results)
    crit(j)=NaN;
    if ~isempty(results(j).theta)
        switch(cr)
            case 2
                crit(j)=results(j).aic;
            case 3
                crit(j)=results(j).bic;
            otherwise
                crit(j) = results(j).deviance;
        end
    end
end

[~,k]=min(crit);
th_ss = results(k).theta;   

figure(1);
present_est_results(y,dt,time,th_ss);



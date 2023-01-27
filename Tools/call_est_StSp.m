function result = call_est_StSp;
% script for estimating a state space model for analyze_GUI.
%
% AUTHOR: dbauer, 31.8.2020.

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
gg = findobj(gcf,'tag','Order_stat');
n = str2num(get(gg,'String'));

gg = findobj(gcf,'tag','crit_stat');
cr =get(gg,'Value');

gg = findobj(gcf,'tag','PEM_stat');
Pbull = 1-get(gg,'Value');



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

result = StSp_I0([y,dt],s,n,floor(sqrt(T*s)),Pbull);
present_est_results(y,dt,time,result.theta);

% compare to existing result -> replace? 
gg= findobj('tag','results_stat');
gres = get(gg,'userdata');

if isempty(gres) % no stationary system estimated up to now -> add. 
    gres = null_model(y,dt);
    gres.time =time;
    gres(n+1)=result;
    set(gg,'userdata',gres);
    return
else % there are results available. 
    if (length(gres)<n+1)
        gres(n+1)=result;
        set(gg,'userdata',gres);
        return;
    else 
        res_cur = gres(n+1); 
        if isempty(res_cur) % system of order n has been estimated. 
            gres(n+1)=result;
            set(gg,'userdata',gres);
            return
        else 
            switch cr 
                case 1
                    crit_new = result.deviance;
                    crit_cur = res_cur.deviance;
                case 2
                    crit_new = result.aic;
                    crit_cur = res_cur.aic;
                case 3
                    crit_new = result.bic;
                    crit_cur = res_cur.bic;
            end
    
            answer = questdlg(sprintf('Use new estimate? Criterion: Current %2.2f. New: %2.2f.',crit_cur,crit_new), 'New estimate','Current','New','Current');
            % Handle response
            switch answer
                case 'New'
                    gres(n+1)=result;
                    set(gg,'userdata',gres);
                    disp('Using new estimate.')
                otherwise 
                    disp('Keeping current system.')
            end
        end
    end
end





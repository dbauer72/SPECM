function result = call_est_SPECM;
% script for estimating a state space model for analyze_GUI.
%
% AUTHOR: dbauer, 31.8.2020.

g = get(gcf,'userdata');
y = g{1};
dt = g{2};
S = g{3};
time = g{4};

[T,s]= size(y);

sdet = 0;
try
    oJoh_j = findobj(gcf,'tag','SJoh_j');
    Joh_j = get(oJoh_j,'value');
    odet = findobj(gcf,'tag','det');
    sdet = get(odet,'value');
catch
    disp('Wrong entry in field detrend.');
    return;
end

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

% get specification for selection
gg = findobj(gcf,'tag','Order');
n = str2num(get(gg,'String'));

gg = findobj(gcf,'tag','crit');
cr =get(gg,'Value');

gg = findobj(gcf,'tag','PEM');
Pbull = 1-get(gg,'Value');

gg = findobj(gcf,'tag','URS_SS');
urs = get(get(gg,'SelectedObject'),'String');

gg = findobj(gcf,'tag','URS_ind_SS');
urs_ind = str2num(get(gg,'String'));

% check for the existence of a stationary estimate with the same order. 
gg_stat = findobj('tag','results_stat');
res_stat = get(gg_stat,'userdata');

restrict.det_res = 0;


% differentiate according to setting of class
switch urs
    case 'I(1)'
        
        if isempty(urs_ind)
            urs_ind = 0;
        end;
        if length(urs_ind)>1
            urs_ind = urs_ind(1);
        end;
%         switch Joh_j
%             case {1,2}
%                 dt = [dt,ones(T,1),[1:T]'];
%             case {3,4}
%                 dt = [dt,1:T]'
%             otherwise
%                 dt = dt; %zeros(T,0);
%         end;
        if (length(res_stat)<n+1)||(~isempty(res_stat(n+1).theta))
            result = SPECM_I1([y,dt],s,n,urs_ind,floor(sqrt(T*s)),Pbull,restrict,res_stat(n+1).theta);
        else
            result = SPECM_I1([y,dt],s,n,urs_ind,floor(sqrt(T*s)),Pbull);
        end
        
    case 'MFI(1)'
        [T,s]= size(y);
        if isempty(urs_ind)
            urs_ind = [0,0,1];
        end;

%         switch Joh_j
%             case {1,2,3,4}
%                 vv = var(dt);
%                 if (min(vv)>10^(-6)) % no constant contained so far -> add. 
%                     dt = [dt,ones(T,1)];
%                 end;
            %otherwise
            %    dt = [zeros(T,0);
%        end;
                
        nmin = sum((urs_ind(:,2)+1).*urs_ind(:,3));
        if (n< nmin)
            warning('SPECM_MFI1: Order specified too small -> changing it.')
            n = nmin;
        end
        %tily = y - dt*(dt\y);
        if (length(res_stat)<n+1)||(~isempty(res_stat(n+1).theta))
            result = SPECM_MFI1([y,dt],s,S,n,urs_ind,floor(sqrt(T*s)),Pbull,restrict,res_stat(n+1).theta);
        else
            result = SPECM_MFI1([y,dt],s,S,n,urs_ind,floor(sqrt(T*s)),Pbull);   
        end
                         
        %result.theta.D = (dt\y)';

    otherwise % I(2) is selected
        [T,s]= size(y);
        if isempty(urs_ind)
            urs_ind = [1,0];
        end;
        
%         switch Joh_j
%             case {1,2,3,4}
%                 const = 1;
%             otherwise
%                 const = 0;
%         end;

        if (length(res_stat)<n+1)||(~isempty(res_stat(n+1).theta))
             result = SPECM_I2([y],s,n,urs_ind,floor(sqrt(T*s)),Pbull,restrict,res_stat(n+1).theta);
        else
             result = SPECM_I2([y],s,n,urs_ind,floor(sqrt(T*s)),Pbull);   
        end       
end;


% compare to existing result -> replace? 
gg= findobj('tag','results');
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
        if isempty(res_cur.theta) % system of order n has been estimated. 
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




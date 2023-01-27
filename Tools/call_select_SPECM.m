function th_ss = call_select_SPECM;
% script for selecting the lag length in the VAR case for analyze_GUI.
%
% OUTPUT:   th_ss ... theta structure of estimate.
%
% AUTHOR: dbauer, 18.8.2020.

g = get(gcf,'userdata');
y = g{1};
s = size(y,2); 
dt = g{2};
S = g{3};
time = g{4};

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

% get specification for selection
gg = findobj(gcf,'tag','MaxO');
nmax = str2num(get(gg,'String'));

gg = findobj(gcf,'tag','crit');
cr =get(gg,'Value');

gg = findobj(gcf,'tag','PEM');
Pbull = 1-get(gg,'Value');



gg = findobj(gcf,'tag','URS_SS');
urs = get(get(gg,'SelectedObject'),'String');

gg = findobj(gcf,'tag','URS_ind_SS');
urs_ind = str2num(get(gg,'String'));

% find out, if a stationary system has been estimated already. 
gg_stat = findobj('tag','results_stat');
res_stat = get(gg_stat,'userdata');

restrict.det_res = 0;

% switch according to selected model 
switch urs
    case 'I(1)' % we want I(1) models to be estimated
        % include deterministics 
        [T,s]= size(y);
        if isempty(urs_ind)
            urs_ind = 0;
        end;
        if length(urs_ind)>1
            urs_ind = urs_ind(1);
        end;
        
        if sdet>0 
            switch Joh_j 
                case {1,2} 
                    dt = [dt,ones(T,1),[1:T]'];
                case {3,4}
                    dt = [dt,[1:T]'];
               % otherwise
               %     dt = zeros(T,0);
            end
        end
        
        results(1) = null_model(y,dt);
        results(1).time =time;

        h=  waitbar(0,'Iterating over order ...','Name','Order estimation');
        for n=1:nmax
            waitbar(n/nmax,h,sprintf('Estimating order %d of %d.',n,nmax));
            if ~isempty(res_stat(n+1))
                [results(n+1)] = SPECM_I1([y,dt],s,n,min(n,urs_ind),floor(sqrt(T)),Pbull,restrict,res_stat(n+1).theta);
            else
                [results(n+1)] = SPECM_I1([y,dt],s,n,min(n,urs_ind),floor(sqrt(T)),Pbull);
            end
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
                par = th2param(result.theta,urs_ind,1);
                s2 = s*(s+1)/2;
                result.param = [results(n).param(1:s2);par;result.theta.D(:)];
                results(n+1) = random_improve_call(result,.01,5);
            end
        end
        close(h);
    case 'MFI(1)' % we want MFI(1) models to be estimated
        [T,s]= size(y);
        if isempty(urs_ind)
            urs_ind = [0,0,1];
        end;

        if sdet>0 
            switch Joh_j 
                case {1,2} 
                    dt = [dt,ones(T,1),[1:T]'];
                case {3,4}
                    dt = [dt,[1:T]'];
               % otherwise
               %     dt = zeros(T,0);
            end
        end
        
        gg= findobj('tag','results');
        gres = get(gg,'userdata');
        
        nmin = sum((urs_ind(:,2)+1).*urs_ind(:,3));
        
        results(1) = null_model(y,dt);
        results(1).time =time;

        h=  waitbar(0,'Iterating over order ...','Name','Order estimation');
        for n=nmin:nmax
            waitbar((n-nmin)/(nmax-nmin),h,sprintf('Estimating order %d of %d.',n,nmax));
            if ~isempty(res_stat(n+1))
                [results(n+1)] = SPECM_MFI1(y,s,S,n,urs_ind,floor(sqrt(T*s)),Pbull,restrict,res_stat(n+1).theta);
            else
                [results(n+1)] = SPECM_MFI1(y,s,S,n,urs_ind,floor(sqrt(T*s)),Pbull);
            end
            results(n+1).time =time;
        end
        close(h);
        
    otherwise % we want I(2) models to be estimated. 
      [T,s]= size(y);
        if isempty(urs_ind)
            urs_ind = [1,0];
        end;
        if length(urs_ind)~=2
            urs_ind = [1,0];
        end;
        if sdet>0 
            switch Joh_j 
                case {1,2} 
                    dt = [dt,ones(T,1),[1:T]'];
                case {3,4}
                    dt = [dt,[1:T]'];
               % otherwise
               %     dt = zeros(T,0);
            end
        end
        
        results(1) = null_model(y,dt);
        results(2) = null_model(y,dt);
        results(1).time =time;
        results(2).time =time;

        c = urs_ind(1)*2+urs_ind(2);
        h=  waitbar(0,'Iterating over order ...','Name','Order estimation');
        % order 1 is left out -> not possible with one I(2) component.
        for n=2:nmax
            if n<c
                results(n)=results(1);
            else
                waitbar(n/nmax,h,sprintf('Estimating order %d of %d.',n,nmax));
                if ~isempty(res_stat(n+1))
                    [results(n+1)] =  SPECM_I2([y],s,n,urs_ind,floor(sqrt(T*s)),Pbull,restrict,res_stat(n+1).theta);
                else
                    [results(n+1)] =  SPECM_I2([y],s,n,urs_ind,floor(sqrt(T*s)),Pbull);
                end
            end
        end
        close(h);  
end


gg= findobj('tag','results');
set(gg,'userdata',results);

pres_select('non-stat');

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



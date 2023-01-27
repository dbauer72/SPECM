function th_est = call_est_VAR;
% script for selecting the lag length in the VAR case for analyze_GUI.
%
% AUTHOR: dbauer, 18.8.2020.

g = get(gcf,'userdata');
y = g{1};
T = size(y,1); 
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

% check, which model to estimate
gg = findobj(gcf,'tag','URS');
urs = get(get(gg,'SelectedObject'),'String');

gg = findobj(gcf,'tag','URS_ind');
urs_ind = str2num(get(gg,'String'));

% get lag length
gg = findobj(gcf,'tag','Lag');
p = str2num(get(gg,'String'));

% levels VAR or VECM?
gg = findobj(gcf,'tag','vecm');
vecm = get(gg,'value');

% now start the engine!
if vecm>0 % VECM estimation 
    % differentiate according to setting of class 
    switch urs
        case 'I(1)'
            gg = findobj(gcf,'tag','Joh_j');
            Joh_j = get(gg,'value');
            [~,~,~,~,~,th_est] = VECM_I1(y,p,size(y,2)-urs_ind(1),Joh_j);
            switch Joh_j 
                case {1,2}
                    u = [ones(T,1),[1:T]'];
                case {3,4} 
                    u = ones(T,1);
                otherwise 
                    u = zeros(T,0);
            end;
        case 'MFI(1)'
            gg = findobj(gcf,'tag','Joh_j');
            Joh_j = get(gg,'value');
            gg = findobj(gcf,'tag','S');
            S = get(gg,'userdata');
            
            [LL,th_est,res,LR] = VECM_MFI1(y,p,S,urs_ind,Joh_j);            
            u=zeros(T,0); 
            
            
        otherwise % I(2) is selected 
            gg = findobj(gcf,'tag','Joh_j');
            Joh_j = get(gg,'value');
            
            [th_est,alpha,beta,zeta,eta,psi,B,Gamma,Omega,et] = VECM_2SI2(y,p,urs_ind,Joh_j);
            if Joh_j<5
                u = ones(T,1);
            else
                u = zeros(T,1);
            end;
    end;
    
else    
    T = size(y,1);
    fdt = dt;
    if ve0>0
        fdt = [fdt,ones(T,1)];
    end
    if ve1>0
        fdt = [fdt,[1:T]'];
    end
    
    if ((ves>0)&&(S>1))
        if ve0>0
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
    
    % protect again double usage of determinstics. 
    [q,r]= qr(fdt);
    dr = diag(r);
    cc = sum(dr>10^(-12));
    fdt = q(:,1:cc);

    % estimate
    u = fdt;
    th_est = VAR(y,u,p);       
    
end

% present estimation results. 
present_est_results(y,u,time,th_est,0,labs);



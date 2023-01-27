function [LRs,rs_est,cv] = call_spec_SPECM;
% script for selecting the order in the VARMA case for analyze_GUI.
%
% AUTHOR: dbauer, 28.8.2020.

g = get(gcf,'userdata');
y = g{1};
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
gg = findobj(gcf,'tag','Order');
n = str2num(get(gg,'String'));

gg = findobj(gcf,'tag','crit');
cr =get(gg,'Value');

gg = findobj(gcf,'tag','PEM');
Pbull = get(gg,'Value');

gg = findobj(gcf,'tag','URS_SS');
urs = get(get(gg,'SelectedObject'),'String');

gg = findobj(gcf,'tag','URS_ind_SS');
urs_ind = str2num(get(gg,'String'));


% add deterministics to dt
[T,s]= size(y);
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
%                 dt = [ones(T,1),[1:T]'];
%             case {3,4}
%                 dt = [1:T]'
%             otherwise
%                 dt = zeros(T,0);
%         end;
        
        gg= findobj('tag','results');
        gres = get(gg,'userdata');
        
        if isempty(gres)
            result = SPECM_I1([y,dt],s,n,urs_ind,floor(sqrt(T*s)),Pbull);
        else
            result = gres(n+1); 
        end        

        
        [LLy,LLn,rs_est,alphahat,betahat,Chat] = test_Johansen_RH_I1(result,0.05,Joh_j,1);
        cv = [];
        % perform estimation for all variants of Joh_j to compare
        % results.
        disp('  Johansen trace test for SPECM form and different determinstics ');
        for j=1:5
            LRs(j,:) = test_Johansen_RH_I1(result,0.05,j,0);
            disp(sprintf(['Det. No. ' num2str(j),': ' repmat('%3.4f ',1,size(LRs,2))],LRs(j,:)));
        end;
        disp('\n\n');
        
    case 'MFI(1)'
%         [T,s]= size(y);
        if isempty(urs_ind)
            urs_ind = [0,0,1];
        end;

        Joh_j = 5;
        
        gg= findobj('tag','results');
        gres = get(gg,'userdata');
        
        if isempty(gres)
            [result] = SPECM_MFI1(y,s,S,n,urs_ind,floor(sqrt(T*s)),Pbull);
        else
            result = gres(n+1); 
        end
        
        thc = result.theta;
        Abar = thc.A-thc.K*thc.C;
                
        [LL] = RH_SPECM_MFI1(y,S,Abar,thc.K,Joh_j);
        [LRs,rs_est] = test_Johansen_SPECM_MFI1(y,Abar,thc.K,0.05,S,Joh_j,1);
        cv = [];
     otherwise % I(2) is selected  
%         [T,s]= size(y);
        if isempty(urs_ind)
            urs_ind = [1,0];
        end;
        if length(urs_ind)~=2
            urs_ind = [1,0];
        end;
         switch Joh_j
             case {1,2,3,4}
                 const = 1;
             otherwise
                 const=0;
         end;
         
        gg= findobj('tag','results');
        gres = get(gg,'userdata');
        
        if isempty(gres)
            result = SPECM_I2([y],s,n,urs_ind,floor(sqrt(T*s)),Pbull);
        else
            result = gres(n+1); 
        end        

        Ae = result.theta.A - result.theta.K*result.theta.C;
        Be = result.theta.K;
        
        % guard against non-minimum-phase solutions. 
        while max(abs(eig(Ae)))>1
            Ae = Ae/2;
        end;
        
        [th_init,Chat,Q1,Q2,alpha,beta,zeta,eta,psi,Best,Gamma] = RH_VECM_LM(y,Ae,Be,urs_ind,const);
        [LRs,rs_est,cv] = test_Johansen_RH_I2(Q1,Q2,s,2000,0.05,0);
        
        % perform estimation for all variants of Joh_j to compare
        % results.
        %disp('  Johansen trace test for SPECM form in I(2) case ');
        LRs
        LRs./cv
        
        figure;
        cc = LRs./cv;
        image((cc>1)*256)

        xlabel('rs(2)');
        ylabel('rs(1)');
        set(gca,'xtick',1:s,'xticklabel',0:(s-1),'ytick',1:s,'yticklabel',0:(s-1));
        title(sprintf('Estimated rs(1): %d, rs(2): %D ',rs_est(1),rs_est(2)));



end;

%present_est_results(y,u,time,result.theta);
 


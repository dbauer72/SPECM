function [th_est,b_eta] = call_spec_VAR;
% script for selecting the lag length in the VAR case for analyze_GUI.
%
% AUTHOR: dbauer, 18.8.2020.

g = get(gcf,'userdata');
y = g{1};
dt = g{2};
S = g{3};
time = g{4};
[T,s] = size(y); 

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
            [LL,b_eta,alphahat,betahat,th_est] = test_Johansen_VECM_I1(y,p,0.05,Joh_j,1); 
                        
            switch Joh_j
                case {1,2} 
                    u = [ones(T,1),[1:T]'];
                case {3,4} 
                    u = ones(T,1);
                otherwise
                    u=zeros(T,0); 
            end; 
            
            % perform estimation for all variants of Joh_j to compare
            % results.
            disp('  Johansen trace test for VECM form and different determinstics ');
            for j=1:5
                LLj(j,:) = test_Johansen_VECM_I1(y,p,0.05,j,0);
                disp(sprintf(['Det. No. ' num2str(j),': ' repmat('%3.4f ',1,size(LLj,2))],LLj(j,:)));
            end;
            disp('');
            
            
        case 'MFI(1)'
            gg = findobj(gcf,'tag','Joh_j');
            Joh_j = get(gg,'value');
            gg = findobj(gcf,'tag','S');
            S = get(gg,'userdata');
            om = [0:(S/2)]'/S*2;
          
            
            [LL,b_eta,th_est,res,LR,ured] = test_Johansen_VECM_MFI1(y,p,0.05,om,S,Joh_j,1);
            Tu = size(ured,1);
            u = [zeros(T-Tu,size(ured,2));ured];

        otherwise % I(2) is selected 
            gg = findobj(gcf,'tag','Joh_j');
            Joh_j = get(gg,'value');
            
            if Joh_j< 5
                const = 1;
            else
                const = 0;
            end;

            [LRs,b_eta,th_est] = test_Johansen_VECM_I2(y,p,0.05,const);
            cv = crit_val_sel(size(y,2),0.05,const);
            
            fprintf('I(2) diagnostics: rows: rank of Pi; cols: rank of Gamma. \n');
            fprintf('Test results/critical values. Values higher than 1 mark rejections.\n')
            LRs./cv
            
            if const
                u = ones(T,1);
            else
                u = zeros(T,1);
            end;

    end;
else % if it is not VECM, the urs does not play a role!

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

    [q,r]= qr(fdt);
    dr = diag(r);
    cc = sum(dr>10^(-12));
    fdt = q(:,1:cc);

    th_est = VAR(y,fdt,p,'LS');
    u=fdt;
    A= th_est.a;
    s = size(y,2);
    switch urs 
        case 'I(1)'
            Pi = eye(s);
            for j=1:p
                Pi = Pi + A(:,j*s+[1:s]);
            end
            figure(3);
            [~,d,~] = svd(Pi);
            plot(sort(abs(diag(d)))*T,'x');
            title('Eigenvalues of estimated matrix Pi*T')

            b_eta = sum(abs(diag(d))>log(T)/T); % very crude estimate of number of common trends!
            fprintf('Crude estimate of the cointegration rank at frequency 0: %d\n',b_eta);
            fprintf('Matrix Pi');
            Pi
        case 'MFI(1)'
            for j=1:S
                fr = 2*(j-1)/S;
                zk = exp(sqrt(-1)*pi*fr);
                Pi = eye(s);
                for j=1:p
                    Pi = Pi + A(:,j*s+[1:s])*zk^j;
                end
                figure(3);
                hold on;
                [~,d,~] = svd(Pi);
                plot(sort(abs(diag(d)))*T,'x');
                title('Eigenvalues of estimated matrix Pi*T')

                b_eta = sum(abs(diag(d))>log(T)/T); % very crude estimate of number of common trends!
                fprintf('Crude estimate of the cointegration rank at frequency %1.2f pi: %d\n',fr,b_eta);
                fprintf('Matrix Pi');
                Pi
            end
        otherwise % I(2) case
            Pi = eye(s);
            for j=1:p
                Pi = Pi + A(:,j*s+[1:s]);
            end
            figure(3);
            [~,d,~] = svd(Pi);
            plot(sort(abs(diag(d)))*T,'x');
            title('Eigenvalues of estimated matrix Pi*T')

            b_eta = sum(abs(diag(d))>log(T)/T); % very crude estimate of number of common trends!
            fprintf('Crude estimate of the cointegration rank at frequency 0: %d\n',b_eta);
            fprintf('Matrix Pi');
            Pi

    end
end

present_est_results(y,u,time,th_est);





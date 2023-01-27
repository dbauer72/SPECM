function [res] = present_est_results(y,u,time,th_est,uncond,labs);
% presents the fit of estimate th_est for dependent variable y and
% deterministics u.
%
% SYNTAX: present_est_results(y,u,th_est);
%
% INPUT: y ... Txs real; observations of y.
%        u ... Txm real; deterministics.
%        th_est ... theta structure of estimates.
%
% OUTPUT: on screen.
%
% AUTHOR: dbauer, 18.8.2020.

% calculate residuals
% two cases: a) pure VARX and conditional likelihood thought -> condition
% omn first p observations.
% b) all else -> convert to state space and run Kalman filter. Starting
% conditions are chosen as zero for all unstalbe directions, acc. to
% stationary distribution else.

if nargin<6
    p = size(y,2);
    labs = cell(p,1);
    for j=1:p
        labs{j} = sprintf('y_%d',j);
    end;
end
if nargin<5
    ee = findobj(gcf,'tag','cond');
    uncond = get(ee,'value');
end

b = th_est.b;
[T,s]=size(y);
if size(b,2)==s
    % this is the pure VARX case.
    a = th_est.a;
    p = round(size(a,2)/s)-1;
    
    % produce matrix of lags
    Z = u(p+1:end,:);
    for j=1:p
        Z = [Z,y(p-j+[1:(T-p)],:)];
    end;
    
    d= th_est.d;
    
    % calculate residuals
    res = y(p+1:end,:) + Z*[-d,a(:,(s+1):end)]';
    Omega = res'*res/(T-p);
    res = [NaN(p,s);res];
    LL = -log(det(Omega))*(T-p)/2 - (T-p)*s/2;
    if (uncond >0)
        LLcon = LL;
        % go to companion form
        Pi = [-a(:,(s+1):end);eye((p-1)*s),zeros(s*(p-1),s)];
        Ci = Pi(1:s,:);
        Bi = [eye(s);zeros(s*(p-1),s)];
        
        % find eigenvalues of large magnitude.
        [t,d] = seig(Pi);
        Pit = d;
        Cit = Ci*t;
        Bit = inv(t)*Bi;
        dd = diag(d);
        if max(abs(dd))<0.99
            Pt = (Bit*Omega*Bit')./(1-dd*dd');
            P = real(t*Pt*t');
            ypt = y([p:-1:1],:)';
            ypv = ypt(:);
            LL = LL - log(det(P))/2 - 0.5*(ypv'*inv(P)*ypv);
        else
            ii = sum(abs(dd)>0.99);
            disp(['Found ',num2str(ii),' unstable components. Assuming zero initial values for these!']);
            Pt = (Bit(ii+1:end,:)*Omega*Bit(ii+1:end,:)')./(1-dd(ii+1:end)*dd(ii+1:end)');
            Ptf = [zeros(ii,p*s);zeros(p*s-ii,ii),Pt]; % covariance of initial state for t=0.
            Q = Bit*Omega*Bit';
            for j=1:p % perform p iterations to get to state at time p.
                Ptf = Pit*Ptf*Pit'+Q;
            end
            Ptf = real(t*Ptf*t');
            ypt = y([p:-1:1],:)';
            ypv = ypt(:);
            
            LL = LL - log(det(Ptf))/2 - 0.5*(ypv'*inv(Ptf)*ypv);
        end
    end
    AIC = -2*LL + 2*(p*s+size(d,2)*s);
    BIC = -2*LL + log(T-p)*(p*s+size(d,2)*s);
else
    if th_est.which == 'SS' % state space form
        %gg = findobj(gcf,'tag','URS');
        %urs = get(get(gg,'SelectedObject'),'String');
        urs = th_est.ur;
        switch urs
            case 'I(1)'
                n = size(th_est.A,1);
                m = size(u,2);
                if ~isempty(th_est.urs)
                    c = th_est.urs(1);
                else
                    c=0;
                end;
                [LL,res] = cal_quasi_like_theta_I1(th_est,[y u],s,m,n,c,0);
                if uncond>0
                    LLcon = LL;
                    [LL,~] = cal_quasi_like_theta_I1(th_est,[y u],s,m,n,c,1);
                end
            case 'MFI(1)'
                n = size(th_est.A,1);
                m = size(u,2);
                urs = th_est.urs;
                [LL,res] = cal_quasi_like_theta_MFI1(th_est,[y u],s,m,n,urs,0);
                if uncond>0
                    LLcon = LL;
                    [LL,~] = cal_quasi_like_theta_I1(th_est,[y u],s,m,n,urs,1);
                end
            case 'I(2)'
                n = size(th_est.A,1);
                m = size(u,2);
                c = th_est.urs;
                
                [LL,res] = cal_quasi_like_theta_I2(th_est,[y u],s,m,n,c,0);
            otherwise % must be stationary then 
                n = size(th_est.A,1);
                m = size(u,2);
                c=0;
                [LL,res] = cal_quasi_like_theta_I1(th_est,[y u],s,m,n,c,0);
                if uncond>0
                    LLcon = LL;
                    [LL,~] = cal_quasi_like_theta_I1(th_est,[y u],s,m,n,c,1);
                end

        end
        Omega = th_est.Omega;
        LL = -LL/2;
    end
end

% put information on screen:

% plot fit

% color scheme
col = {'k','b','r','m','c','g','y'}; 

figure(4);
hold off;
ax1 = subplot(1,2,1);
hold on;
labs2 = cell(0,1);

for j=1:size(y,2)
    plot(time,y(:,j),[col{1+rem(j-1,7)}]);
    plot(time,y(:,j)-res(:,j),[col{1+rem(j-1,7)},'-.']);
    labs2{end+1} = labs{j};
    labs2{end+1} = [labs{j},':pred'];    
end;
legend(labs2); 
title('Obs. and predictions');

ax2 = subplot(1,2,2);
hold off
for j=1:size(y,2)
    plot(time,res(:,j),[col{1+rem(j-1,7)}]);
    hold on;
end
legend(labs); 

title('Residuals');

linkaxes([ax1 ax2],'x')

% print basic facts.
fprintf('Sample size: %d.  Dimension: %d\n',T,s);
fprintf('Likelihood: %6.2f. Deviance: %6.2f\n',LL,-2*LL);
AIC = -2*LL+2*th_est.num_param;
BIC = -2*LL+log(T)*th_est.num_param;
fprintf('AIC:  %2.6f. BIC: %2.6f\n',AIC,BIC);
if uncond>0
    fprintf('Cond. Likelihood: %6.2f. Cond. Deviance: %6.2f\n',LLcon,-2*LLcon);
end;

fprintf('Omega (from final estimates)')
Omega
fprintf('Omega (from selection)')
th_est.Omega


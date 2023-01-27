function fore_plot_call(fig,Test,ini,os); 
% fore_plot_call executes the calculations and plots for fore_plot.
%
% SYNTAX: fore_plot_call(fig,Test,ini,os); 
%
% INPUT:  fig ... handle to figure
%         Test ... real; point to start prediction
%         ini  ... binary; 1: estimate effects due to initial values.
%         os   ... binary; 1: one step predictions, else starting from
%         Test+1 on. 
%
% OUTPUT: plots.
%
% REMARK: data is contained in est_results vector in userdata of figure.
% 
% AUTHOR: dbauer, 3.11.2020

figure(fig);
results = get(fig,'userdata');


col = {'k','b','r','m','c','g','y'}; 

hold off;
ax1 = subplot(1,2,1);
hold on;

title('Obs. and predictions');

ax2 = subplot(1,2,2);
hold off

title('Residuals');

linkaxes([ax1 ax2],'x')


Test = max(1,round(Test));
Tval = Test;

% pick the model wanted
obj = findobj(gcf,'tag','Nr');
ch = get(obj,'value');

% obtain prediction horizon
obj = findobj(gcf,'tag','ph');
ph_st = get(obj,'String');

try 
    ph = round(str2num(ph_st));
catch 
    ph = 1;
end

ph = max(1,ph);
ph = min(ph,Tval);

%for j=1:length(results)
res = results(ch);
y = res.y;

[T,s] = size(y);
dt = res.dt;
m = size(dt,2);
th = res.theta;
n = size(th.A,1);
if ~isprop(res,'time')
    time = 1:T;
else
    time = res.time;
end
if isempty(time)
    time = 1:T;
end

% extract initial values, of ini
%     if ini>0
%         OT = th.C;
%         for t=2:T
%             OT = [th.C;OT*th.A];
%         end
%         xx = [OT,kron(dt,eye(s))];
%         vy = reshape(y',s*T,1);
%        tvy = vy - xx*(xx(1:(Test*s),:)\vy(1:(s*Test),:));
%        ty = reshape(tvy,s,T)';
%     else
if isempty(th.D)
    th.D = zeros(s,m);
    th.B = zeros(n,m);
end

ty = y - dt*th.D';
%     end
% calculate Kalman filter on estimation part

% initial variance
ur = res.ur;
urs = res.urs;
switch ur
    case 'I(1)'
        c = urs;
    case 'I(2)'
        c = sum(urs);
    case 'MFI(1)'
        c = sum((urs(:,2)+1).*urs(:,3));
    otherwise
        c = 0;
end

n = res.n;
nb = n-c;
Abull = th.A(c+1:n,c+1:n);
C = th.C;
K = th.K;
A = th.A;
Omega = th.Omega;
Q = K*Omega*K';
Qb = Q(c+1:n,c+1:n);
pb = inv(eye(nb^2)-kron(Abull',Abull'))*Qb(:);
Pbull = reshape(pb,nb,nb);

P0 = [0*eye(c),zeros(c,n-c);zeros(n-c,c),Pbull];

P0 = (P0+P0')/2; % P(1|0).
% --- initialize the Kalman filter ---
x0= zeros(n,1); % x(1|0)

% --- run the Kalman filter ---
tres = ty*0;
hres = tres;
tres(1,:)=ty(1,:); % e(1)
hres(1:ph,:) = tres(1:ph,:); % for ph step ahead prediction, first ph residuals are the output.
Omegat= C*P0*C'+Omega; % Omega(1|0)
Kt = (A*P0*C'+K*Omegat)*inv(Omegat);
xf = A*x0 + Kt*tres(1,:)';
Pkg1 = A*P0*A' + Q- Kt*Omegat*Kt'; % P(2|1)


if os
    Test = T;
end
for t=2:(Test) % filter equations
    Omegat = C*Pkg1*C'+Omega;
    iOm = inv(Omegat);
    tres(t,:)= ty(t,:) - xf'*C';
    Kt = (A*Pkg1*C'+K*Omega)*iOm;
    xf = A*xf+ Kt*tres(t,:)';

    % generate ph-step ahead prediction. 
    if (t+ph<=T)
        hres(t+ph-1,:)= ty(t+ph-1,:)-xf'*(A')^(ph-1)*C';    
    end
    Pkg1 = (A * Pkg1 *A') +Q - (Kt*Omegat*Kt');
    Pkg1 = (Pkg1+Pkg1')/2;
end

% split into one step prediction or full prediction on validation data.
if os == 0 % full prediction for rest of sample
    % continue to predict
    xf = A^(ph-1)*xf;
    for t=Test+1:(T-ph)
        hres(t,:)= ty(t+ph-1,:) - xf'*C';
        xf = A*xf;
    end
end

% produce plots
axes(ax1);
hold off
S = res.S;
for m=1:s
    plot(time,y(:,m),col{m});
    hold on;
    plot(time,y(:,m)-hres(:,m),[col{m},'-.']);
end
plot(time(Tval)*[1,1],[min(y(:)),max(y(:))],'k');
set(ax1,'xlim',[time(1),time(T)]);
title('Obs. and predictions');
axes(ax2);
hold off;
for m=1:s
    plot(time,hres(:,m),col{m});
    hold on;
end
set(ax2,'xlim',[time(1),time(T)]);
title('Residuals');
% write out results
rmse_e = sqrt(sum(sum(hres(1:Tval,:).^2))/Tval);
rmse_v = sqrt(sum(sum(hres(Tval+1:T,:).^2))/(T-Tval));

fprintf('Model: %d: RMSE in sample: %f; RMSE out of sample: %f\n\n',j,rmse_e,rmse_v);
%end




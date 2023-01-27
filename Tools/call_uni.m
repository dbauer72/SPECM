function ok = call_uni(j);
% callback for univariate analyses, support of analyze GUI.
%
% SYNTAX:  call_uni(j);
%   
% INPUT:  j ... integer; coordinate of series to analyze. 
%
% AUTHOR: dbauer, 17.8.2020
g = get(gcf,'userdata');
y = g{1};
dt = g{2};
S = g{3};
time = g{4};

% 
d = findobj(gcf,'tag',sprintf('D%d',j));
try 
    dd = str2num(get(d,'string')); 
catch
    disp('Wrong entry in field diff.');
    return;
end

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


lab = get(findobj(gcf,'tag',sprintf('L%d',j)),'string');
ycur = y(:,j);
T0 = length(ycur);

% seasonal differences? 
ds = 0;
if S>1 
   tds = findobj(gcf,'tag','DS');
   ds = get(tds,'value');
end

for l=1:dd
    if ds
        ycur = ycur((S+1):end)-ycur(1:(end-S));
        time = time((S+1):end);
        lab = sprintf('SD(%s)',lab);
    else
        ycur = diff(ycur);
        lab = sprintf('D(%s)',lab);
        time = time(2:end);
    end
    
end

% detrend variables 
T = length(ycur);
if ~isempty(dt)
    fdt = dt(end-length(ycur)+1:end,:);
else
    fdt = zeros(length(ycur),0);
end

if ve0>0 
    fdt = [fdt,ones(T,1)];
end
if ve1>0 
    fdt = [fdt,[1:T]'];
end

if ((ves>0)&&(S>1))
    for j=2:S
        dums = zeros(T,1);
        dums([j:S:T])=1;
        fdt = [fdt,dums];
    end
end

% orthogonalize in order to prevent 'doubling up'
[q,r]= qr(fdt);
dr = abs(diag(r));
cc = sum(dr>10^(-12));
fdt = q(:,1:cc);

yhat = fdt*(fdt\ycur);
ycur = ycur- yhat;

% plot the data set
% identify, if plot already exists 
obj = findobj('MenuBar','figure');
strs = {obj(:).Name};
ind = strmatch('P/ACF',strs);
if ~isempty(ind)
    figure(obj(ind(1)).Number);
else
    figure('Name','P/ACF')
end
pacf_sample(ycur,10,lab,1,yhat,time);

% output on adf test
maxL = max(round(sqrt(T)),round(T/10));
nlag = max(1,aicest(ycur,1,maxL)); % always include at least one lagged difference in adf test.

% choose deterministics
p=-1;
if ve1>0
    p=1;
else
    if ve0>0
        p=0;
    end
end

if (ds == 0)||(S ~= 4)
    % adf test
    results = adf(ycur,p,nlag);
    prt(results);
    adf_v = results.adf;
    cr = results.crit;
    if cr(2)>adf_v
        disp('ADF test rejects existence of unit root at 5%%!');
    else
        disp('ADF test does not reject at 5%% -> variable might be integrated.');
    end
else
%     % HEGY test for seasonal integration for quarterly observations.
%     %            det = 1, no deterministic components
% %                = 2, constant (default)
% %                = 3, constant & 3 seasonal dummies
% %                = 4, constant & trend
% %                = 5, constant & 3 seasonal dummies & trend
%     p=5;    
%     if ve0+ve1+ves == 0
%         p = 1;
%     end
%     if (ve0==1)&&(ve1+ves==0)
%         p=2;
%     end
%     if (ve1== 0)&&(ve0==1)&&(ves==1)
%         p=3;
%     end
%     if (ve1+ve0==2)&&(ves==0)
%         p=4;
%     end;
%     % HEGY constrained to sample size smaller than 200 -> cut off for
%     % longer time series
%     hegy(ycur(1:min(200,length(ycur))),0.05,p,nlag);
    % alternative: Portmanteau test. 
    [tP,pp,z] = test_season_port(ycur,1+floor(S/2),2*nlag); % use two times the lag length of an autoregression as M. 
    for jj=1:length(tP)
        if pp(jj)>0.05
            if (z(jj,2)==1)
                fprintf('Portmanteau test does not reject (p-value: %2.2f) \n the existence of unit root at frequency %2.2f (unit root at z=%2.2f+i*%2.2f) at 5%%! \n',pp(jj),z(jj,1),real(exp(sqrt(-1)*z(jj,1))),imag(exp(sqrt(-1)*z(jj,1))));
            else 
                fprintf('Portmanteau test does not reject (p-value: %2.2f) \n the existence of unit root at frequency %2.2f (unit root at z=%2.2f) at 5%%! \n',pp(jj),z(jj,1),real(exp(sqrt(-1)*z(jj,1))));
            end
        else
            if (z(jj,2)==1)
                fprintf('Portmanteau test rejects (p-value: %2.2f) \n the existence of unit root at frequency %2.2f (unit root at z=%2.2f+i*%2.2f) at 5%%! \n',pp(jj),z(jj,1),real(exp(sqrt(-1)*z(jj,1))),imag(exp(sqrt(-1)*z(jj,1))));
            else 
                fprintf('Portmanteau test rejects (p-value: %2.2f) \n the existence of unit root at frequency %2.2f (unit root at z=%2.2f) at 5%%! \n',pp(jj),z(jj,1),real(exp(sqrt(-1)*z(jj,1))));
            end
        end
    end
end





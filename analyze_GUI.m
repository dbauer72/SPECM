function ok = analyze_GUI(y,labs,S,dt,time);
% analyze_GUI provides access to a number of analyses using a GUI.
%
% SYNTAX: analyze_GUI(y);
%
% INPUT: y ... T x s real matrix of observations.
%        labs ... s x 1 cell vector of names for observation coordinates.
%        S ... integer; number of observations per year. 
%
% OUTPUT: none, only graphical. 
%
% REMARK: Caution, outputs many analyses to figures and the workspace.
%
% AUTHOR: dbauer, 3.4.2020.


if nargin<3
    S=1
end;

figure;

[T,s]= size(y);
if nargin<5
    time = [1:T]/S;
end;


if nargin<4
    dt = zeros(T,0);
end

if nargin<2
    labs = cell(s,1);
    for j=1:s
        labs{j} = sprintf('y%d',j);
    end;
end;

set(gcf,'userdata',{y,dt,S,time,labs});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% univariate analyses 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uni = uipanel('Title','Univariate','FontSize',12,'BackgroundColor','white','Position',[.01 .8 .98 .18]);

uicontrol('Parent',uni,'style','text','String',sprintf('S: %d',S),'userdata',S,'tag','S','units','normalized','Position',[(0)/(s+2),.7,0.95/(s+2),.25]);

for j=1:s
    uicontrol('Parent',uni,'style','text','String',labs{j},'tag',sprintf('L%d',j),'units','normalized','Position',[(j)/(s+2),.7,1/(s+2),.25]);
    uicontrol('Parent',uni,'style','edit','String','0','tag',sprintf('D%d',j),'units','normalized','Position',[(j)/(s+2),.4,1/(s+2),.25]);
    uicontrol('Parent',uni,'style','pushbutton','String','Analyze!','units','normalized','Position',[(j)/(s+2),.05,1/(s+2),.25],'call',sprintf('call_uni(%d)',j));    
end;

uicontrol('Parent',uni,'style','text','String','Diff.','units','normalized','Position',[0.05,.4,1/(s+2)-0.05,.25]);

% detrending 
uicontrol('Parent',uni,'style','text','String','Const.','units','normalized','Position',[(s+1)/(s+2),.7,0.8/(s+2),.25]);
uicontrol('Parent',uni,'style','text','String','lin.tr.','units','normalized','Position',[(s+1)/(s+2),.4,0.8/(s+2),.25]);
uicontrol('Parent',uni,'style','text','String','seas.dumm.','units','normalized','Position',[(s+1)/(s+2),.05,0.8/(s+2),.25]);

uicontrol('Parent',uni,'style','checkbox','value',1,'tag','T0','units','normalized','Position',[(s+1.8)/(s+2),.7,0.2/(s+2),.20]);
uicontrol('Parent',uni,'style','checkbox','value',1,'tag','T1','units','normalized','Position',[(s+1.8)/(s+2),.4,0.2/(s+2),.20]);
uicontrol('Parent',uni,'style','checkbox','value',1,'tag','TS','units','normalized','Position',[(s+1.8)/(s+2),.05,0.2/(s+2),.20]);

% if S>1: introduce seasonal differencing. 
if S>1
    uicontrol('Parent',uni,'style','text','String','Seas.Diff.','units','normalized','Position',[0.05,.05,1/(s+2)-0.05,.25]);
    uicontrol('Parent',uni,'style','checkbox','value',0,'tag','DS','units','normalized','Position',[1/(s+2)-0.05,.05,0.05,.25]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multivariate VAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%

var = uipanel('Title','VAR','FontSize',12,'BackgroundColor','white','Position',[.01 .6 .98 .18]);

s=8; 
uicontrol('Parent',var,'style','text','String','Maxlag','units','normalized','Position',[1/(s+2),.6,0.9/(s+2),.35]);
uicontrol('Parent',var,'style','edit','String','10','tag','MaxL','units','normalized','Position',[2/(s+2),.6,0.9/(s+2),.35]);
uicontrol('Parent',var,'style','text','String','Lag','units','normalized','Position',[1/(s+2),.1,0.9/(s+2),.35]);
uicontrol('Parent',var,'style','edit','String','1','tag','Lag','units','normalized','Position',[2/(s+2),.1,0.9/(s+2),.35]);

sel ='th = call_select_VAR;';
uicontrol('Parent',var,'style','pushbutton','String','Select lag!','units','normalized','Position',[3/(s+2),.1,0.9/(s+2),.8],'call',sel);    

est = 'th = call_est_VAR;';
uicontrol('Parent',var,'style','pushbutton','String','Estimate!','units','normalized','Position',[4/(s+2),.1,0.9/(s+2),.8],'call',est);    

uicontrol('Parent',var,'style','text','String','Uncond. Like','units','normalized','Position',[5/(s+2),.6,0.9/(s+2),.35]);
uicontrol('Parent',var,'style','checkbox','value',0,'tag','cond','units','normalized','Position',[5/(s+2),.1,0.9/(s+2),.35]);

% VECM specs
uicontrol('Parent',var,'style','text','String','VECM','units','normalized','Position',[6/(s+2),.7,0.9/(s+2),.25]);
uicontrol('Parent',var,'style','checkbox','value',0,'tag','vecm','units','normalized','Position',[6/(s+2),.4,0.9/(s+2),.25]);
uicontrol('Parent',var,'style','popupmenu','value',5,'String',{'unr. c+t','rest. c+t','unr. c','rest. c','none'},'tag','Joh_j','units','normalized','Position',[6/(s+2),.1,0.9/(s+2),.25]);


bg = uibuttongroup(var,'Position',[7/(s+2),.1,0.9/(s+2),.8],'tag','URS');
              
% Create three radio buttons in the button group.
uicontrol(bg,'Style','radiobutton','String','I(1)','units','normalized','Position',[.1 .7 .8 .25]);
uicontrol(bg,'Style','radiobutton','String','MFI(1)','units','normalized','Position',[.1 .4 .8 .25]);
uicontrol(bg,'Style','radiobutton','String','I(2)','units','normalized','Position',[.1 .1 .8 .25]);

uicontrol('Parent',var,'style','text','String','URS','units','normalized','Position',[8/(s+2),.7,0.9/(s+2),.25]);
uicontrol('Parent',var,'style','edit','String','[1]','tag','URS_ind','units','normalized','Position',[8/(s+2),.4,0.9/(s+2),.25]);

spec = '[th_est,b_eta] = call_spec_VAR;';
uicontrol('Parent',var,'style','pushbutton','String','Specify!','units','normalized','Position',[8/(s+2),.1,0.9/(s+2),.25],'call',spec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stationary State Space
% modeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mod_st = uipanel('Title','Stationary State Space Model','FontSize',12,'BackgroundColor','white','Position',[.01 .4 .98 .18]);

% similar setup as with VARs
s=8; 
ggn =uicontrol('Parent',mod_st,'style','Pushbutton','String','Pres. Sel.','tag','results_stat','units','normalized','Position',[0/(s+2),.6,0.9/(s+2),.35],'call','pres_select(''stat'');');
results_stat(1) = null_model(y,dt);
results_stat(1).time =time;
set(ggn,'userdata',results_stat);

% replace entry

use_result_stat = ['gg= findobj(''tag'',''results_stat'');results = get(gg,''userdata'');gn = findobj(gcf,''tag'',''Order_stat'');n = str2num(get(gn,''String''));',...
    'if isempty(results),clear results;end;results(n+1)=result;set(gg,''userdata'',results);'];

uicontrol('Parent',mod_st,'style','Pushbutton','String','Use result','units','normalized','Position',[0/(s+2),.1,0.9/(s+2),.35],'call',use_result_stat);
uicontrol('Parent',mod_st,'style','text','String','Max order','units','normalized','Position',[1/(s+2),.6,0.9/(s+2),.35]);
uicontrol('Parent',mod_st,'style','edit','String','10','tag','MaxO_stat','units','normalized','Position',[2/(s+2),.6,0.9/(s+2),.35]);
uicontrol('Parent',mod_st,'style','text','String','Order','units','normalized','Position',[1/(s+2),.1,0.9/(s+2),.35]);
uicontrol('Parent',mod_st,'style','edit','String','1','tag','Order_stat','units','normalized','Position',[2/(s+2),.1,0.9/(s+2),.35]);

sel_stat ='results = call_select_StSp;';
uicontrol('Parent',mod_st,'style','pushbutton','String','Select order!','units','normalized','Position',[3/(s+2),.3,0.9/(s+2),.6],'call',sel_stat);  
uicontrol('Parent',mod_st,'style','popupmenu','value',1,'String',{'Deviance','AIC','BIC'},'tag','crit_stat','units','normalized','Position',[3/(s+2),.05,0.9/(s+2),.2]);


est_stat = 'result = call_est_StSp;';
imp_stat = 'result = call_improve_StSp;';

uicontrol('Parent',mod_st,'style','pushbutton','String','Improve!','units','normalized','Position',[4/(s+2),.6,0.9/(s+2),.35],'call',imp_stat);    

uicontrol('Parent',mod_st,'style','pushbutton','String','Estimate!','units','normalized','Position',[4/(s+2),.1,0.9/(s+2),.35],'call',est_stat); 

uicontrol('Parent',mod_st,'style','text','String','PEM','units','normalized','Position',[5/(s+2),.6,0.9/(s+2),.35]);
uicontrol('Parent',mod_st,'style','checkbox','value',0,'tag','PEM_stat','units','normalized','Position',[5/(s+2),.1,0.9/(s+2),.35]);

uicontrol('Parent',mod_st,'style','text','String','Det. Input','units','normalized','Position',[6/(s+2),.7,0.9/(s+2),.25]);
uicontrol('Parent',mod_st,'style','checkbox','value',0,'tag','det_stat','units','normalized','Position',[6/(s+2),.4,0.9/(s+2),.25]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-stationary State Space
% modeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



mod = uipanel('Title','Non-stationary State Space Model','FontSize',12,'BackgroundColor','white','Position',[.01 .2 .98 .18]);

% similar setup as with VARs
s=8; 
ggn =uicontrol('Parent',mod,'style','Pushbutton','String','Pres. Sel.','tag','results','units','normalized','Position',[0/(s+2),.6,0.9/(s+2),.35],'call','pres_select(''non-stat'');');
results(1) = null_model(y,dt);
results(1).time =time;
set(ggn,'userdata',results);

% replace entry

use_result = ['gg= findobj(''tag'',''results'');results = get(gg,''userdata'');gn = findobj(gcf,''tag'',''Order'');n = str2num(get(gn,''String''));',...
    'if isempty(results),clear results;end;results(n+1)=result;set(gg,''userdata'',results);'];

uicontrol('Parent',mod,'style','Pushbutton','String','Use result','units','normalized','Position',[0/(s+2),.1,0.9/(s+2),.35],'call',use_result);
uicontrol('Parent',mod,'style','text','String','Max order','units','normalized','Position',[1/(s+2),.6,0.9/(s+2),.35]);
uicontrol('Parent',mod,'style','edit','String','10','tag','MaxO','units','normalized','Position',[2/(s+2),.6,0.9/(s+2),.35]);
uicontrol('Parent',mod,'style','text','String','Order','units','normalized','Position',[1/(s+2),.1,0.9/(s+2),.35]);
uicontrol('Parent',mod,'style','edit','String','1','tag','Order','units','normalized','Position',[2/(s+2),.1,0.9/(s+2),.35]);

sel ='results = call_select_SPECM;';
uicontrol('Parent',mod,'style','pushbutton','String','Select order!','units','normalized','Position',[3/(s+2),.3,0.9/(s+2),.6],'call',sel);  
uicontrol('Parent',mod,'style','popupmenu','value',1,'String',{'Deviance','AIC','BIC'},'tag','crit','units','normalized','Position',[3/(s+2),.05,0.9/(s+2),.2]);


est = 'result = call_est_SPECM;';
imp = 'result = call_improve_SPECM;';

uicontrol('Parent',mod,'style','pushbutton','String','Improve!','units','normalized','Position',[4/(s+2),.6,0.9/(s+2),.35],'call',imp);    

uicontrol('Parent',mod,'style','pushbutton','String','Estimate!','units','normalized','Position',[4/(s+2),.1,0.9/(s+2),.35],'call',est); 

uicontrol('Parent',mod,'style','text','String','PEM','units','normalized','Position',[5/(s+2),.6,0.9/(s+2),.35]);
uicontrol('Parent',mod,'style','checkbox','value',0,'tag','PEM','units','normalized','Position',[5/(s+2),.1,0.9/(s+2),.35]);

uicontrol('Parent',mod,'style','text','String','Det. Input','units','normalized','Position',[6/(s+2),.7,0.9/(s+2),.25]);
uicontrol('Parent',mod,'style','checkbox','value',0,'tag','det','units','normalized','Position',[6/(s+2),.4,0.9/(s+2),.25]);
uicontrol('Parent',mod,'style','popupmenu','value',5,'String',{'unr. c+t','rest. c+t','unr. c','rest. c','none'},'tag','SJoh_j','units','normalized','Position',[6/(s+2),.1,0.9/(s+2),.25]);

% Create three radio buttons in the button group.
bg = uibuttongroup(mod,'Position',[7/(s+2),.1,0.9/(s+2),.8],'tag','URS_SS');             

uicontrol(bg,'Style','radiobutton','String','I(1)','units','normalized','Position',[.1 .7 .8 .25]);
uicontrol(bg,'Style','radiobutton','String','MFI(1)','units','normalized','Position',[.1 .4 .8 .25]);
uicontrol(bg,'Style','radiobutton','String','I(2)','units','normalized','Position',[.1 .1 .8 .25]);

uicontrol('Parent',mod,'style','text','String','URS','units','normalized','Position',[8/(s+2),.7,0.9/(s+2),.25]);
uicontrol('Parent',mod,'style','edit','String','[1]','tag','URS_ind_SS','units','normalized','Position',[8/(s+2),.4,0.9/(s+2),.25]);

spec = '[th_est,b_eta] = call_spec_SPECM;';
uicontrol('Parent',mod,'style','pushbutton','String','Specify!','units','normalized','Position',[8/(s+2),.1,0.9/(s+2),.25],'call',spec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimated models 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thetas = uipanel('Title','Inspect models','FontSize',12,'BackgroundColor','white','Position',[.01 .01 .98 .18]);

% what do we want to inspect? 

% impulse responses 
plot_impulse = 'M = str2num(get(findobj(''tag'',''imp_lags''),''String''));impulse_plot(ths,M);';
uicontrol('Parent',thetas,'style','pushbutton','String','Impulse','units','normalized','Position',[0.05,.6,0.2,.35],'call',plot_impulse);
uicontrol('Parent',thetas,'style','edit','String','10','tag','imp_lags','units','normalized','Position',[0.05,.05,0.2,.35]);

% poles / zeros
plot_pz = 'polezero_plot(ths);';
uicontrol('Parent',thetas,'style','pushbutton','String','Pole/Zero','units','normalized','Position',[0.3,.6,0.2,.35],'call',plot_pz);

% forecasts 
plot_fore = 'fore_plot(result);';
uicontrol('Parent',thetas,'style','pushbutton','String','Forecast','units','normalized','Position',[0.55,.6,0.2,.35],'call',plot_fore);
uicontrol('Parent',thetas,'style','edit','String','4','tag','fore_lags','units','normalized','Position',[0.55,.05,0.2,.35]);

% hypothesis tests. 
test_hyp = 'hypothesis_GUI(result);';
uicontrol('Parent',thetas,'style','pushbutton','String','Hypothesis Testing on "result"','units','normalized','Position',[0.8,.05,0.15,.9],'call',test_hyp);

function ok = fore_plot(resu);
% fore_plot plots forecaasts for all models contained in resu
%
% SYNTAX: fore_plot(resu);
%
% INPUT: resu ... vector of est_result objects.
%    
% OUTPUT: plot (GUI).
%
% AUTHOR: dbauer, 3.11.2020

figure; 
y = resu(1).y;
T=size(y,1);
Test = T-20;

set(gcf,'userdata',resu);

%  setting up GUI.
plot_call = 'ini = get(findobj(''tag'',''Ini''),''value'');os = get(findobj(''tag'',''OS''),''value'');Tval = get(findobj(''tag'',''Tval''),''value'');fore_plot_call(gcf,Tval,ini,os);';
% text
uicontrol('style','text','string','Results refer to model #','units','normalized','position',[0.05,0.95,0.15,0.05]);
uicontrol('style','text','string','Est. Init. Eff.?','units','normalized','position',[0.05,0.0,0.1,0.05]);
uicontrol('style','text','string','Rolling pred.?','units','normalized','position',[0.2,0.0,0.1,0.05]);
uicontrol('style','text','string','Pred. horizon?','units','normalized','position',[0.35,0.0,0.1,0.05]);
uicontrol('style','text','string','Split Est. data-Val. data','units','normalized','position',[0.55,0.0,0.1,0.05]);
%
%% controls
numbers = cell(0);
for j=1:length(resu)
    numbers{j}= num2str(j);
end;

uicontrol('style','popupmenu','String',numbers,'tag','Nr','units','normalized','position',[0.2,0.95,0.05,0.05],'value',1,'call',plot_call);
uicontrol('style','checkbox','tag','Ini','units','normalized','position',[0.15,0.0,0.05,0.05],'value',1);
uicontrol('style','checkbox','tag','OS','units','normalized','position',[0.3,0.0,0.05,0.05],'value',1);
uicontrol('style','edit','tag','ph','units','normalized','position',[0.4,0.0,0.05,0.05],'String',num2str(1));
uicontrol('style','Slider','tag','Tval','units','normalized','position',[0.65,0.0,0.3,0.05],'min',1,'max',T,'sliderstep',[.01,.1],'value',Test,'call',plot_call);

fore_plot_call(gcf,Test,1,1); 
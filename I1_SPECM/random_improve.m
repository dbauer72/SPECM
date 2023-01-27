function ok = random_improve(result);
%
%


figure(17);

set(gcf,'Position',[258   430   272   149]);

set(gcf,'userdata',{result});

% call string
r_i_call = ['g = get(gcf,''userdata'');result=g{1};ll_cur = result.deviance; ',... 
    'k_obj = findobj(gcf,''tag'',''k'');k = str2num(get(k_obj,''string''));',...
    'M_obj = findobj(gcf,''tag'',''M'');M = fix(str2num(get(M_obj,''string'')));',...
    'result = random_improve_call(result,k,M);g{1}=result;set(gcf,''userdata'',g);',...
    'ol_obj = findobj(gcf,''tag'',''ov_v'');set(ol_obj,''string'',sprintf(''%4.2f'',ll_cur));',...
    'nl_obj = findobj(gcf,''tag'',''nv_v'');set(nl_obj,''string'',sprintf(''%4.2f'',result.deviance));']; 

%  setting up GUI.
uicontrol('style','edit','string','0.01','tag','k','units','normalized','position',[0.05,0.05,0.3,0.3]);
uicontrol('style','edit','string','5','tag','M','units','normalized','position',[0.55,0.05,0.3,0.3]);
uicontrol('style','pushbutton','string','Run optimization','units','normalized','position',[0.05,0.8,0.9,0.2],'call',r_i_call);

% info on improvement
uicontrol('style','text','string','Old val.','tag','ov','units','normalized','position',[0.05,0.6,0.3,0.2]);
uicontrol('style','text','string','New val.','tag','nv','units','normalized','position',[0.55,0.6,0.3,0.2]);
uicontrol('style','text','string','NA','tag','ov_v','units','normalized','position',[0.05,0.35,0.4,0.2]);
uicontrol('style','text','string','NA','tag','nv_v','units','normalized','position',[0.55,0.35,0.4,0.2]);



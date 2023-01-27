function ok = hypothesis_GUI(result)
% hypothesis_GUI provides access to a number of hypothesis tests using a GUI.
%
% SYNTAX: hypothesis_GUI(result);
%
% INPUT: result ... results structure from estimation. 
%
% OUTPUT: none, only graphical. 
%
% AUTHOR: dbauer, 16.2.2022.

if ~strmatch(class(result),'est_result') 
    message('Input arument must be an est_result object.');
    return
end

figure(17);
set(gcf,'userdata',result);

switch result.theta.ur
    case 'I(1)' % I(1) case 
       ANZ = 12; 
       % text fields 
       uni = uipanel('Title','I(1)','FontSize',12,'BackgroundColor','white','Position',[.01 .05 .98 .95]);

       uicontrol('Parent',uni,'style','text','String','Deterministics:','units','normalized','Position',[.7,1-1/(ANZ),.1,1/ANZ]);

       uicontrol('Parent',uni,'style','text','String',sprintf('s: %d, r: %d',result.s,result.s-result.urs),'units','normalized','Position',[.05,1-1/(ANZ),.25,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','Hypotheses on beta','units','normalized','Position',[.05,1-2/(ANZ),.25,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','beta = beta_0','units','normalized','Position',[.15,1-3/(ANZ),.3,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','beta = H phi','units','normalized','Position',[.15,1-4/(ANZ),.3,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','beta = [b,b_perp varphi]','units','normalized','Position',[.15,1-5/(ANZ),.3,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','Hypotheses on alpha','units','normalized','Position',[.05,1-6/(ANZ),.25,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','type A: alpha = A psi','units','normalized','Position',[.15,1-7/(ANZ),.3,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','type B: alpha =[a,a_perp psi]','units','normalized','Position',[.15,1-8/(ANZ),.3,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','Hypotheses on both alpha and beta','units','normalized','Position',[.05,1-9/(ANZ),.25,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','type C: alpha = A psi, beta = H phi','units','normalized','Position',[.15,1-10/(ANZ),.3,1/ANZ]);

       uicontrol('Parent',uni,'style','text','String','\beta_0:','units','normalized','Position',[.5,1-3/(ANZ),.15,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','H:','units','normalized','Position',[.5,1-4/(ANZ),.15,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','b:','units','normalized','Position',[.5,1-5/(ANZ),.15,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','A:','units','normalized','Position',[.5,1-7/(ANZ),.15,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','a:','units','normalized','Position',[.5,1-8/(ANZ),.15,1/ANZ]);

       % edit fields 
       uicontrol('Parent',uni,'style','popupmenu','value',5,'String',{'unr. c+t','rest. c+t','unr. c','rest. c','none'},'tag','HJoh_j','units','normalized','Position',[.8,1-1/(ANZ),.15,1/ANZ]);
       
       uicontrol('Parent',uni,'style','edit','tag','beta','String','[]','units','normalized','Position',[.7,1-3/(ANZ),.3,1/ANZ]);
       uicontrol('Parent',uni,'style','edit','tag','H','String','[]','units','normalized','Position',[.7,1-4/(ANZ),.3,1/ANZ]);
       uicontrol('Parent',uni,'style','edit','tag','b','String','[]','units','normalized','Position',[.7,1-5/(ANZ),.3,1/ANZ]);
       uicontrol('Parent',uni,'style','edit','tag','A','String','[]','units','normalized','Position',[.7,1-7/(ANZ),.3,1/ANZ]);
       uicontrol('Parent',uni,'style','edit','tag','a','String','[]','units','normalized','Position',[.7,1-8/(ANZ),.3,1/ANZ]);

       % radio button groups. 
       bg1 = uibuttongroup('Parent',uni,'Visible','off','units','normalized','Position',[.05,1-5/(ANZ),.1,3/ANZ],'tag','type_beta');
              
        % Create three radio buttons in the button group.
        r1 = uicontrol(bg1,'Style','radiobutton','String','A','units','normalized','Position',[0,.7,1,.3]);
        r2 = uicontrol(bg1,'Style','radiobutton','String','B','units','normalized','Position',[0,.4,1,.3]);
        r3 = uicontrol(bg1,'Style','radiobutton','String','C','units','normalized','Position',[0,.1,1,.3]);
        bg1.Visible = 'on';

       % radio button groups. 
       bg2 = uibuttongroup('Parent',uni,'Visible','off','units','normalized','Position',[.05,1-8/(ANZ),.1,2/ANZ],'tag','type_alpha');
              
        % Create three radio buttons in the button group.
        r12 = uicontrol(bg2,'Style','radiobutton','String','A','units','normalized','Position',[0,.5,1,.5]);
        r22 = uicontrol(bg2,'Style','radiobutton','String','B','units','normalized','Position',[0,0,1,.5]);
              
        bg2.Visible = 'on';

        % call backs
        test_beta = ['result = get(gcf,''userdata'');y = result.y;theta= result.theta; Abar = theta.A-theta.K*theta.C;rest = ''y'';' ,...
            'obj = findobj(''tag'',''type_beta''); butt = get(obj,''SelectedObject''); restrict.type =get(butt,''String'');' ,...
            'obj = findobj(''tag'',''beta''); restrict.beta = str2num(get(obj,''String''));' ,...
            'obj = findobj(''tag'',''H''); restrict.H = str2num(get(obj,''String''));' ,...
            'obj = findobj(''tag'',''b''); restrict.b = str2num(get(obj,''String''));' ,...
            'obj = findobj(''tag'',''HJoh_j''); Joh_j = get(obj,''value'');' ,...
            '[LL,~,~,df] = RH_specm_H0(y,Abar,theta.K,result.s-result.urs,rest,Joh_j,restrict); '];


        test_alpha = ['result = get(gcf,''userdata'');y = result.y;theta= result.theta; Abar = theta.A-theta.K*theta.C;rest = ''y'';' ,...
            'obj = findobj(''tag'',''type_alpha''); butt = get(obj,''SelectedObject''); restrict.type =get(butt,''String'');' ,...
            'obj = findobj(''tag'',''A''); restrict.A = str2num(get(obj,''String''));' ,...
            'obj = findobj(''tag'',''a''); restrict.a = str2num(get(obj,''String''));' ,...
            'obj = findobj(''tag'',''HJoh_j''); Joh_j = get(obj,''value'');' ,...
            '[LL,~,~,df] = RH_specm_Ha(y,Abar,theta.K,result.s-result.urs,rest,Joh_j,restrict); '];
 
        test_alphabeta = ['result = get(gcf,''userdata'');y = result.y;theta= result.theta; Abar = theta.A-theta.K*theta.C;rest = ''y'';' ,...
            'restrict.type =''C'';' ,...
            'obj = findobj(''tag'',''A''); restrict.A = str2num(get(obj,''String''));' ,...
            'obj = findobj(''tag'',''H''); restrict.H = str2num(get(obj,''String''));' ,...
            'obj = findobj(''tag'',''HJoh_j''); Joh_j = get(obj,''value'');' ,...
            '[LL,~,~,df] = RH_specm_Ha(y,Abar,theta.K,result.s-result.urs,rest,Joh_j,restrict); '];

        % action buttons 
        uicontrol('Parent',uni,'style','Pushbutton','String','Test beta!','units','normalized','Position',[.7,1-2/ANZ,.2,1/ANZ],'call',test_beta);
        uicontrol('Parent',uni,'style','Pushbutton','String','Test alpha!','units','normalized','Position',[.7,1-6/ANZ,.2,1/ANZ],'call',test_alpha);
        uicontrol('Parent',uni,'style','Pushbutton','String','Test alpha and beta!','units','normalized','Position',[.7,1-9/ANZ,.2,1/ANZ],'call',test_alphabeta);
        
    case 'MFI(1)' % MFI(1) case 
       uni = uipanel('Title','MFI(1): Complex Roots','FontSize',12,'BackgroundColor','white','Position',[.01 .7 .98 .3]);

       ANZ=4;
       uicontrol('Parent',uni,'style','text','String','Deterministics:','units','normalized','Position',[.7,1-1/(ANZ),.1,1/ANZ]);
       uicontrol('Parent',uni,'style','popupmenu','value',5,'String',{'unr. c+t','rest. c+t','unr. c','rest. c','none'},'tag','HJoh_j','units','normalized','Position',[.8,1-1/(ANZ),.15,1/ANZ]);
       
       uicontrol('Parent',uni,'style','text','String','Real valuedness hypotheses for beta_j.','units','normalized','Position',[.01,1-1/(ANZ),.5,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','Frequency','units','normalized','Position',[.01,1-2/(ANZ),.18,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','Real valuedness hypotheses for alpha_j.','units','normalized','Position',[.01,1-3/(ANZ),.5,1/ANZ]);
       uicontrol('Parent',uni,'style','text','String','Frequency','units','normalized','Position',[.01,1-4/(ANZ),.18,1/ANZ]);

       urs = result.urs; 
       ind = find(urs(:,2)==1); % find all complex unit root frequencies.
       ni = length(ind);
       for j=1:length(ind)
            uicontrol('Parent',uni,'style','text','String',sprintf('%2.3f',urs(ind(j),1)),'units','normalized','Position',[.2+0.8*2*(j-1)/(2*ni),1-2/(ANZ),1/(2*ni),1/ANZ]);
            uicontrol('Parent',uni,'style','checkbox','value',0,'units','normalized','Position',[.2+0.8*(2*(j-1)+1)/(2*ni),1-2/(ANZ),1/(2*ni),1/ANZ],'tag','B','userdata',ind(j));
       end
        for j=1:length(ind)
            uicontrol('Parent',uni,'style','text','String',sprintf('%2.3f',urs(ind(j),1)),'units','normalized','Position',[.2+0.8*2*(j-1)/(2*ni),1-4/(ANZ),1/(2*ni),1/ANZ]);
            uicontrol('Parent',uni,'style','checkbox','value',0,'units','normalized','Position',[.2+0.8*(2*(j-1)+1)/(2*ni),1-4/(ANZ),1/(2*ni),1/ANZ],'tag','A','userdata',ind(j));
       end

       uicontrol('Parent',uni,'style','text','String','Active: ','units','normalized','Position',[.55,1-1/ANZ,.15,1/ANZ]);
       uicontrol('Parent',uni,'style','checkbox','value',0,'units','normalized','Position',[.7,1-1/ANZ,.05,1/ANZ],'tag','act_com');

       % call back 

       hStems = findobj('-regexp','Tag','^(?!Steady State$).');

       test_com = ['result = get(gcf,''userdata'');y = result.y;theta= result.theta; Abar = theta.A-theta.K*theta.C;',...
           'objsA = findobj(''-regexp'',''tag'',''^A\w*'',''-and'',''checked'',''on'');objsB = findobj(''-regexp'',''tag'',''^B\w*'',''-and'',''checked'',''on'');',...
           'betas = [];alphas = [];for j=1:length(objsA),alphas = [alphas,(get(objsA(j),''userdata'')];end;',...
           'for j=1:length(objsB),betas = [betas,(get(objsA(j),''userdata'')];end;',...
           'if isempty(betas), restrict.type = ''B'';end, if isempty(alphas), restrict.type = ''A'';end, if ~isempty(alphas)&(~isempty(betas)), restrict.type =''C'';end,',...
           'restrict.betafr = betas; restrict.alphafr = alphas;',...
           'obj = findobj(''tag'',''HJoh_j''); Joh_j = get(obj,''value'');' ,...
           '[LLA,df] = RH_SPECM_MFI1_H(y,result.S,Abar,theta.K,Joh_j,result.urs,restrict);'];

       test_z1 = ['result = get(gcf,''userdata'');y = result.y;theta= result.theta; Abar = theta.A-theta.K*theta.C;' ,...
            'obj = findobj(''tag'',''type_beta_z1''); butt = get(obj,''SelectedObject''); restrict.type =[''I1.1.'',get(butt,''String'')];' ,...
            'obj = findobj(''tag'',''beta0_z1''); restrict.beta = str2num(get(obj,''String''));' ,...
            'obj = findobj(''tag'',''H_z1''); restrict.H = str2num(get(obj,''String''));' ,...
            'obj = findobj(''tag'',''b_z1''); restrict.b = str2num(get(obj,''String''));' ,...
           'obj = findobj(''tag'',''HJoh_j''); Joh_j = get(obj,''value'');' ,...
           '[LLA,df] = RH_SPECM_MFI1_H(y,result.S,Abar,theta.K,Joh_j,result.urs,restrict);'];

       test_zm1 = ['result = get(gcf,''userdata'');y = result.y;theta= result.theta; Abar = theta.A-theta.K*theta.C;' ,...
            'obj = findobj(''tag'',''type_beta_zm1''); butt = get(obj,''SelectedObject''); restrict.type =[''I1.m1.'',get(butt,''String'')];' ,...
            'obj = findobj(''tag'',''beta0_zm1''); restrict.beta = str2num(get(obj,''String''));' ,...
            'obj = findobj(''tag'',''H_zm1''); restrict.H = str2num(get(obj,''String''));' ,...
            'obj = findobj(''tag'',''b_zm1''); restrict.b = str2num(get(obj,''String''));' ,...
            'obj = findobj(''tag'',''HJoh_j''); Joh_j = get(obj,''value'');' ,...
           '[LLA,df] = RH_SPECM_MFI1_H(y,result.S,Abar,theta.K,Joh_j,result.urs,restrict);'];

       test_all = ['obj = findobj(''tag'',''act_com''); if get(obj,''checked'',''on''),',test_com,'end;',...
           'obj = findobj(''tag'',''act_z1''); if get(obj,''checked'',''on''),',test_z1,'end;',...
           'obj = findobj(''tag'',''act_zm1''); if get(obj,''checked'',''on''),',test_zm1,'end;']; 

       uicontrol('Parent',uni,'style','pushbutton','String','Test all!','units','normalized','Position',[.75,1-1/ANZ,.25,1/ANZ],'call',test_all);

        % panel for tests at z=1.
       uni2 = uipanel('Title','MFI(1): z=1','FontSize',12,'BackgroundColor','white','Position',[.01 .05 .48 .65]);
       
       
       ANZ2 = 4;
       uicontrol('Parent',uni2,'style','text','String','beta = beta_0: ','units','normalized','Position',[.31,1-1/(ANZ2),.3,1/ANZ2]);
       uicontrol('Parent',uni2,'style','text','String','beta = H phi. H:','units','normalized','Position',[.31,1-2/(ANZ2),.3,1/ANZ2]);
       uicontrol('Parent',uni2,'style','text','String','beta = [b,b_perp varphi]. b:','units','normalized','Position',[.31,1-3/(ANZ2),.3,1/ANZ2]);
       bg1 = uibuttongroup('Parent',uni2,'Visible','off','units','normalized','Position',[.01,1/ANZ2,.3,3/ANZ2],'tag','type_beta_z1');
              
       uicontrol('Parent',uni2,'style','edit','tag','beta_z1','String','[]','units','normalized','Position',[.7,1-1/(ANZ2),.3,1/ANZ2]);
       uicontrol('Parent',uni2,'style','edit','tag','H_z1','String','[]','units','normalized','Position',[.7,1-2/(ANZ2),.3,1/ANZ2]);
       uicontrol('Parent',uni2,'style','edit','tag','b_z1','String','[]','units','normalized','Position',[.7,1-3/(ANZ2),.3,1/ANZ2]);

        % Create three radio buttons in the button group.
        r1 = uicontrol(bg1,'Style','radiobutton','String','A','units','normalized','Position',[0,.7,1,.3]);
        r2 = uicontrol(bg1,'Style','radiobutton','String','B','units','normalized','Position',[0,.4,1,.3]);
        r3 = uicontrol(bg1,'Style','radiobutton','String','C','units','normalized','Position',[0,.1,1,.3]);
        bg1.Visible = 'on';

       uicontrol('Parent',uni2,'style','text','String','Active: ','units','normalized','Position',[.01,0,.3,1/ANZ2]);
       uicontrol('Parent',uni2,'style','checkbox','value',0,'units','normalized','Position',[.31,0,.3,1/ANZ2],'tag','act_z1');

        % panel for tests at z=-1.
       uni3 = uipanel('Title','MFI(1): z=-1','FontSize',12,'BackgroundColor','white','Position',[.51 .05 .48 .65]);
       uicontrol('Parent',uni3,'style','text','String','beta = beta_0: ','units','normalized','Position',[.31,1-1/(ANZ2),.3,1/ANZ2]);
       uicontrol('Parent',uni3,'style','text','String','beta = H phi. H:','units','normalized','Position',[.31,1-2/(ANZ2),.3,1/ANZ2]);
       uicontrol('Parent',uni3,'style','text','String','beta = [b,b_perp varphi]. b:','units','normalized','Position',[.31,1-3/(ANZ2),.3,1/ANZ2]);
       bg3 = uibuttongroup('Parent',uni3,'Visible','off','units','normalized','Position',[.01,1/ANZ2,.3,3/ANZ2],'tag','type_beta_zm1');
              
       uicontrol('Parent',uni3,'style','edit','tag','beta_zm1','String','[]','units','normalized','Position',[.7,1-1/(ANZ2),.3,1/ANZ2]);
       uicontrol('Parent',uni3,'style','edit','tag','H_zm1','String','[]','units','normalized','Position',[.7,1-2/(ANZ2),.3,1/ANZ2]);
       uicontrol('Parent',uni3,'style','edit','tag','b_zm1','String','[]','units','normalized','Position',[.7,1-3/(ANZ2),.3,1/ANZ2]);

        % Create three radio buttons in the button group.
        r1 = uicontrol(bg3,'Style','radiobutton','String','A','units','normalized','Position',[0,.7,1,.3]);
        r2 = uicontrol(bg3,'Style','radiobutton','String','B','units','normalized','Position',[0,.4,1,.3]);
        r3 = uicontrol(bg3,'Style','radiobutton','String','C','units','normalized','Position',[0,.1,1,.3]);
        bg3.Visible = 'on';

       uicontrol('Parent',uni3,'style','text','String','Active: ','units','normalized','Position',[.01,0,.3,1/ANZ2]);
       uicontrol('Parent',uni3,'style','checkbox','value',0,'units','normalized','Position',[.31,0,.3,1/ANZ2],'tag','act_zm1');


    case 'I(2)' % I(2) case 
        msgbox('No hypothesis test for I(2) processes are implemented as of now.')
    otherwise % must be stationary then 
        msgbox('Hypothesis testing works on unit root processes!'); 
end

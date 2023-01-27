function result = call_improve_SPECM;
% script for improving a state space model for analyze_GUI.
%
% AUTHOR: dbauer, 31.8.2020.

g = get(gcf,'userdata');
y = g{1};
dt = g{2};
S = g{3};
time = g{4};

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

gg= findobj('tag','results');
gres = get(gg,'userdata');

if isempty(gres)
    disp('No results obtained so far, nothing to improve!');
    return
else

    % differentiate according to setting of class
    try
        result = gres(n+1);
        if isempty(result)
            disp('No result available for this order!')
            return
        end
    catch
        disp('No result available for this order!')
        return
    end

    random_improve(result);
end

gres(n+1)=result;

set(gg,'userdata',gres);



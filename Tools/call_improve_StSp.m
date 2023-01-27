function result = call_improve_StSp;
% script for improving a state space model for analyze_GUI.
%
% AUTHOR: dbauer, 18.1.2023.

g = get(gcf,'userdata');
y = g{1};
dt = g{2};
S = g{3};
time = g{4};

try
    odet = findobj(gcf,'tag','det_stat');
    sdet = get(odet,'value');
catch
    disp('Wrong entry in field deterministic.');
    return;
end

% get specification for selection
gg = findobj(gcf,'tag','Order_stat');
n = str2num(get(gg,'String'));

gg = findobj(gcf,'tag','crit_stat');
cr =get(gg,'Value');

gg = findobj(gcf,'tag','PEM_stat');
Pbull = 1-get(gg,'Value');


gg= findobj('tag','results_stat');
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
        result = NaN; 
        disp('No result available for this order!')
        return
    end

    random_improve(result);
end

gres(n+1)=result;

set(gg,'userdata',gres);



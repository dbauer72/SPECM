function ok = impulse_plot(thc,M,labs);
% prepares an interactive  plot of the impulse response function for multivariate transfer functions
% 
% SYNTAX:  impulse_plot(thc,M,labs);
%
% INPUT:  thc ... array of theta structures.
%         M   ... time horizon for impulse response. 
%         labs ... labels for the entries of thv. 
%
% OUTPUT: none.
%
% AUTHOR: dbauer, 14.1.2020.
%

fig = figure; 
fig.Position(3:4) = [440 320];

if nargin<2
    M= 100;
end;

m = length(thc);

% differentiate, if a vector of theta_structures or a vector of results is
% supplied. 

if isa(thc(1),'est_result')
    for j=1:m
        ths(j) = thc(j).theta;
    end
else
    for j=1:m
        ths(j) = poly2ss_th(thc(j));
    end
end
s = size(ths(1).C,1);

% make sure all have the same number of inputs
mmax = 0;
for j=1:m
    mmax = max(mmax,size(ths(j).D,2));
end

for j=1:m
    mm =size(ths(j).D,2);
    if (mm<mmax)
        if (mm>0)
            ths(j).D(:,(mm+1):mmax)=0;
            ths(j).B(:,(mm+1):mmax)=0;        
        else
            ths(j).D = zeros(s,mmax);
            ths(j).B = zeros(size(ths(j).A,1),mmax);
        end
    end
end


for j=1:m
    ir{j} = impulse(ths(j),M);
    if nargin<3
        labs{j} = sprintf('th: %d',j);
    end
end
for j=1:s
    vals_o{j}= num2str(j);
    vals_i{j}= ['I:',num2str(j)];
end;

for j=1:size(ths(1).D,2)
    vals_i{s+j} = ['E:',num2str(j)];    
end

i=1;j=1;
hold('on');
for jj=1:m
    irc = ir{jj};
    p(jj) = plot(0:size(irc,3)-1,squeeze(irc(i,j,:)));
end;

legend(labs);

set(gcf,'userdata',{ir,p});
%  setting up GUI.
plot_call = ['g = get(gcf,''userdata'');ir = g{1};p=g{2};m = length(p);obj1 = findobj(''tag'',''In'');ic = get(obj1,''Value'');obj2 = findobj(''tag'',''Out'');jc = get(obj2,''Value'');for j=1:m,irc = ir{j};bc = squeeze(irc(jc,ic,:));set(p(j),''YData'',bc);end;'];

% text
uicontrol('style','text','string','Channel','units','normalized','position',[0.05,0.9,0.1,0.1]);
uicontrol('style','text','string','Out:','units','normalized','position',[0.25,0.9,0.1,0.1]);
uicontrol('style','text','string','In:','units','normalized','position',[0.55,0.9,0.1,0.1]);
%
%% drop down 
uicontrol('style','popupmenu','string',vals_o,'tag','Out','units','normalized','position',[0.35,0.9,0.1,0.1],'call',plot_call);
uicontrol('style','popupmenu','string',vals_i,'tag','In','units','normalized','position',[0.65,0.9,0.1,0.1],'call',plot_call);
%
%% plot 
%uicontrol('style','pushbutton','string','Plot!','units','normalized','position',[0.85,0.9,0.1,0.1],'call',plot_call);

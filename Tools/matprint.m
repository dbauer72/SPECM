function ok = matprint(M,stdM,ylab,xlab)
% provides printout of matrix with standard deviations in parenthesis. 
%
% SYNTAX: mprint(M,stdM,ylab,xlab);
%
% INPUT: M ... mxn real matrix.
%        stdM ... m*n x 1 real vector of standard deviations.
%        ylab, xlab ... cell arrays of strings indicating the labels for
%        coordinates. 
%
% OUTPUT: on screen.
% 
% AUTHOR: dbauer, 9.3.2022

if nargin< 3
    xlab = {};
    ylab = {};
end
if nargin <4
    xlab = {};
end;

xstr = [];
if ~isempty(xlab)
    nx = length(xlab{1});
    for j=2:length(xlab);
        nx = max(nx,length(xlab{j}));
    end
    nx = max(nx,22);
    for j=1:length(xlab)
        xlab{j}((end+1):nx) = ' ';
        xstr = [xstr,'   ',xlab{j}];
    end;
end

if ~isempty(ylab)
    ny = length(ylab{1});
    for j=2:length(ylab)
        ny = max(ny,length(ylab{j}));
    end

    for j=1:length(ylab)
        ylab{j}((end+1):ny) = ' ';
    end
else
     for j=1:size(M,1)
        ylab{j}= ' ';
    end
end

% reformat stds 
[m,n]=size(M);
try 
    sM = reshape(stdM,m,n);
catch
    error('Wrong dimension of stds.')
end

if ~isempty(xstr)
    fprintf('                %s\n',xstr);
end
% now cycle over rows 
for a=1:m
    str = [ylab{a},'   '];
    for b = 1:n
        str = [str,cprint(10,M(a,b)),' (',cprint(10,sM(a,b)),')  '];
    end
    fprintf([str,'\n']);
end


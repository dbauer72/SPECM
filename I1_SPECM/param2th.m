function [th,par,dth] = param2th(par,s,n,c,restrict);
% converts parameters into th structures via system matrices.
% same as param2syst but without D, Omega.
% 
% SYNTAX: [th,par,dth] = param2th(par,s,n,c,restrict); 
%
% INPUT: param ... dx1 parameter vector.
%        s     ... integer; dimension of endogenous vars.
%        n     ... integer; state dimension.
%        c     ... integer; number of common trends.
%        restrict ... structure describing the restrictions. 
%
% OUTPUT: th ... theta structure corresponding to state space system.
%         par ... remaining parameters for D. 
%         dth      ... if wanted in the output contains the derivative of
%                       the system matrices with respect to the parameters. 
%                      Array: dth.A: array nxnxd of derivatives of A w.r.t.
%                      the parameters. 
%
% REMARK: uses the parameterization corr. to the param paper. 
%
% AUTHOR: dbauer, 27.10.2019.

if nargin<5
    restrict = [];
end

if nargout>2
    np =length(par);
    dA = zeros(n,n,np);
    dK = zeros(n,s,np);
    dC = zeros(s,n,np);
end

% --- find number of params in C1 and K1 ---
if isfield(restrict,'H')
    nparc =restrict.nparc;
    parc = par(1:nparc);
    par(1:nparc)=[];
    C1 = feval(restrict.par2ortho,parc,s,c,restrict.H);
else
    nparc = s*c-c*(c+1)/2;
    parc = par(1:nparc);
    par(1:nparc)=[];
    C1 = par2ortho_LR(parc,s,c);
    if nargout>2 
        for j=1:length(parc)
            dC(:,1:c,j) = dpar2ortho_LR(parc,s,c,j);
        end
    end
end

npark = s*c - c*(c-1)/2;
cur =nparc; 

% --- fill in parameters in K1 ---
park = par(1:npark);
K1 = zeros(c,s);
for j=1:c
    K1(j,j:end)=park(1:(s-j+1));
    park(1:(s-j+1))=[];
    if nargout>2
        for jj=1:(s-j+1)
            dK(j,j+jj-1,cur+jj)= 1;
        end
        cur = cur+s-j+1;
    end
end;

par(1:npark)=[];
cur = nparc+npark;

% --- remaining parameters for stationary part ---
parbull = par(1:(n-c)*s*2);
if nargout>2
    [Abull,Kbull,Cbull,dAKC]= param2mat_bull(parbull,s,n-c);
else
    [Abull,Kbull,Cbull]= param2mat_bull(parbull,s,n-c);
end
par(1:(n-c)*s*2)=[];

% --- fill in the submatrices ---
A = [eye(c),zeros(c,n-c);zeros(n-c,c),Abull];
K = [K1;Kbull];
C = [C1,Cbull];

% --- fill into theta structure ---
th = theta_urs;
th.which = 'SS';
th.A = A;
th.K=K;
th.C=C;

if nargout>2
    % fill in 
    for j=(cur+1):(np-length(par))
        dA(:,:,j) = [zeros(c,n);zeros(n-c,c),squeeze(dAKC.A(:,:,j-cur))];
        dC(:,:,j) = [zeros(s,c),squeeze(dAKC.C(:,:,j-cur))];
        dK(:,:,j)= [zeros(c,s);squeeze(dAKC.K(:,:,j-cur))];
    end  
    
    % write into dth
    dth.A = dA;
    dth.K = dK;
    dth.C = dC; 
end

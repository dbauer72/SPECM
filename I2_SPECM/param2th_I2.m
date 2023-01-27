function [th,par,dth] = param2th_I2(par,s,n,ds,restrict);
% converts parameters into th structures via system matrices.
% same as param2syst_I2 but without D, Omega.
% 
% SYNTAX: [th,par,dth] = param2th_I2(par,s,n,ds,restrict); 
%
% INPUT: param ... dx1 parameter vector.
%        s     ... integer; dimension of endogenous vars.
%        n     ... integer; state dimension.
%        ds     ... pair of integers; number of common trends: ds(1) I(2),
%                       ds(2) additional I(1) trends.  
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
%         par ... [parc,park,parbullet,parD,parOmega].
%
% AUTHOR: dbauer, 6.2.2020.

if nargin<5
    restrict = [];
end

if nargout>2
    np =length(par);
    dA = zeros(n,n,np);
    dK = zeros(n,s,np);
    dC = zeros(s,n,np);
end

c = 2*ds(1)+ds(2);

% --- find number of params in C1, C2, C3  ---

% --- C1 ---
nparc1 = s*ds(1)-ds(1)*(ds(1)+1)/2;
parc1 = par(1:nparc1);
par(1:nparc1)=[];
[C1,QL] = par2ortho_LR(parc1,s,ds(1));
if nargout>2
    for j=1:length(parc1)
        dC(:,1:c,j) = dpar2ortho_LR(parc1,s,c,j);
    end
end

% --- C3: C3'C1 = 0, C3'C3 = I ---
nparc3 = (s-ds(1))*ds(2)-ds(2)*(ds(2)+1)/2;
parc3 = par(1:nparc3);
par(1:nparc3)=[];
[C3h,QL2] = par2ortho_LR(parc3,s-ds(1),ds(2));
C3 = QL*[zeros(ds(1),ds(2));C3h];

% --- C2: C2'[C1,C3]=0 ---
nparc2 = ds(1)*(s-sum(ds));
parc2 = par(1:nparc2);
par(1:nparc2)=[];
Lambda= reshape(parc2,s-sum(ds),ds(1));
C2 = QL*[zeros(ds(1),ds(1));[QL2*[zeros(ds(2),ds(1));Lambda]]];

% --- count params ---
cur =nparc1+nparc2+nparc3; 


% --- K1, K2, K3 ---
npark1 = s*ds(1);
park1 = par(1:npark1);
par(1:npark1)=[];
npark2 = s*ds(1) - ds(1)*(ds(1)-1)/2;
park2 = par(1:npark2);
par(1:npark2)=[];
npark3 = s*ds(2) - ds(2)*(ds(2)-1)/2;
park3 = par(1:npark3);
par(1:npark3)=[];

K1 = reshape(park1,ds(1),s);
K2 = K1*0;
for j=1:ds(1)
    K2(j,j:end)=park2(1:(s-j+1)); 
    park2(1:(s-j+1))=[];
end

K3 = zeros(ds(2),s);
for j=1:ds(2)
    K3(j,j:end)=park3(1:(s-j+1)); 
    park3(1:(s-j+1))=[];
end

Ku = [K1;K2;K3]; 

cur =npark1+npark2+npark3; 


% --- remaining parameters for stationary part ---
parbull = par(1:(n-c)*s*2);
if nargout>2
    [Abull,Kbull,Cbull,dAKC]= param2mat_bull(parbull,s,n-c);
else
    [Abull,Kbull,Cbull]= param2mat_bull(parbull,s,n-c);
end
par(1:(n-c)*s*2)=[];

% --- fill in the submatrices ---
A = [eye(ds(1)),eye(ds(1)),zeros(ds(1),n-2*ds(1));zeros(ds(1)+ds(2),ds(1)),eye(ds(1)+ds(2)),zeros(ds(1)+ds(2),n-c);zeros(n-c,c),Abull];
K = [Ku;Kbull];
C = [C1,C2,C3,Cbull];

% --- fill into theta structure ---
th = theta_urs();
th.which = 'SS';
th.A = A;
th.K=K;
th.C=C;
th.ur = 'I(2)';

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

function [th,par,dth] = param2th_MFI1(par,s,n,urs,restrict);
% converts parameters into th structures via system matrices.
% same as param2syst but without D, Omega.
% 
% SYNTAX: [th,par,dth] = param2th(par,s,S,n,urs,restrict); 
%
% INPUT: param ... dx1 parameter vector.
%        s     ... integer; dimension of endogenous vars.
%        n     ... integer; state dimension.
%        urs   ... unit root structure. 
%        restrict ... structure describing the restrictions. 
%
% OUTPUT: th ... theta structure corresponding to state space system.
%         par ... remaining parameters for D. 
%         dth      ... if wanted in the output contains the derivative of
%                       the system matrices with respect to the parameters. 
%                      Array: dth.A: array nxnxd of derivatives of A w.r.t.
%                      the parameters. 
%
% REMARK: + uses the parameterization corr. to the param paper. 
%         + adds parameterization of joint matrix C_j for all unit roots.  
%
% AUTHOR: dbauer, 14.10.2020

sizOm = s*(s+1)/2;

Omega = fill_lowtri(par(1:sizOm),s);

par(1:sizOm) = [];
if nargin<4
    restrict = [];
end

if nargout>2 % gradient wanted?
    np =length(par);
    dA = zeros(n,n,np);
    dK = zeros(n,s,np);
    dC = zeros(s,n,np);
end

if nargin<5
    restrict.det_res=0;
end
% --- fill in the parameters in Cu ---
C = zeros(s,n); 
K = C';
A = zeros(n,n); 
cur = 0;

% --- add different parameterization for equal C_j for all unit roots ---
equalC =0;
if isfield(restrict,'equal')
    if restrict.equal == 1
        equalC = 1;
        ck = urs(1,3);%ck must be the same then for all k.
        nparc = (s-ck)*ck+ck*(ck-1)/2;
        Ck = par2ortho_LR(par(1:nparc),s,ck); 
        par(1:nparc) = [];
    end
end


for j=1:size(urs,1) % cycle over unit roots 
    ck = urs(j,3);
    zk = exp(sqrt(-1)*urs(j,1)*pi);
    if urs(j,2)==0 % real root 
        % ----- Ck -------
        nparc = (s-ck)*ck+ck*(ck-1)/2;
        if equalC == 0
            Ck = par2ortho_LR(par(1:nparc),s,ck);
            par(1:nparc) = [];
        end 
        C(:,cur+[1:ck])=Ck;
        % ----- Ak -------
        A(cur+[1:ck],cur+[1:ck])=zk*eye(ck);
        
        % ----- Kk -------
        Kk = zeros(ck,s);
        Kk([1:ck],1:ck)=diag(par(1:ck)); % diagonal.
        par(1:ck)=[];
        for a=1:ck % upper diagonal.
            for b=a+1:ck
                Kk(a,b)= par(1);
                par(1)=[];
            end;
        end
        
        np = ck*(s-ck);
        kval = par(1:np);
        par(1:np)=[];
        
        Kk(:,ck+1:end) = reshape(kval,ck,s-ck);
        K(cur+[1:ck],:)=real(Kk);
        
        cur = cur+ck;
    else
        % ---------- Ck ------
        if equalC == 0
            nparc = 2*s*ck-ck^2;
            Ck = par2ortho_LR_c(par(1:nparc),s,ck);
            C(:,cur+[1:2*ck]) = [real(Ck),imag(Ck)];
            par(1:nparc) = [];
        else 
            nparc = ck^2;
            omk = par2ortho_LR_c(par(1:nparc),ck,ck);
            Ckcomp = Ck*omk;
            C(:,cur+[1:2*ck]) = [real(Ckcomp),imag(Ckcomp)];
            par(nparc)=[];
        end;
        % --------- Ak -----
        A(cur+[1:2*ck],cur+[1:2*ck])=kron([real(zk),imag(zk);-imag(zk),real(zk)],eye(ck));
        
        % --------- Kk ------
        Kk = zeros(ck,s);
        Kk([1:ck],1:ck)=diag(par(1:ck)); % diagonal.
        par(1:ck)=[];
        for a=1:ck % upper diagonal.
            for b=a+1:ck
                Kk(a,b)= par(1)+sqrt(-1)*par(2);
                par(1:2)=[];
            end;
        end
        
        % remaining block column. 
        np = ck*(s-ck);
        kval = par(1:np) + sqrt(-1)*par(np+[1:np]);
        par(1:2*np)=[];
        
        Kk(:,ck+1:end) = reshape(kval,ck,s-ck);
        K(cur+[1:ck],:)=real(Kk)*2;
        K(cur+ck+[1:ck],:)=imag(Kk)*2;
        
        cur = cur+2*ck;
    end
    
end
   


% --- remaining parameters for stationary part ---
parbull = par(1:(n-cur)*s*2);
if nargout>2
    [Abull,Kbull,Cbull,dAKC]= param2mat_bull(parbull,s,n-cur);
else
    [Abull,Kbull,Cbull]= param2mat_bull(parbull,s,n-cur);
end
par(1:(n-cur)*s*2)=[];

% --- rest for D ----
m = floor(length(par)/s);
D = reshape(par(1:s*m),s,m);

% --- fill in the submatrices ---
A(cur+1:end,cur+1:end) = Abull;
K(cur+1:end,:) = Kbull;
C(:,cur+1:end) = Cbull;

% --- fill into theta structure ---
th = theta_urs;
th.which = 'SS';
th.A = real(A);
th.K=real(K);
th.C=real(C);
th.Omega = real(Omega); 
th.D = D; 
th.B = zeros(n,m); 

if nargout>3 % gradient wanted. 
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

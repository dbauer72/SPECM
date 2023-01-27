function [Abull,Kbull,Cbull,dAKC]= param2mat_bull(parbull,s,n);
% writes the parameters into the matrices for system in echelon form.
%
% SYNTAX: [Abull,Kbull,Cbull,dAKC]= param2mat_bull(parbull,s,n);
%
% INPUT: parbull ... dx1 parameter vector.
%        s       ... integer; output dimension.
%        n       ... integer; state dimension.
%
% OUTPUT: Abull,Kbull,Cbull ... system matrices.
%         dAKC   ... if nargout>3, then also derivatives of the system
%         matrices with respect to the parameters are calculated. 
%                     dAKC.A: array of dimension nxnxd. 
%                     dAKC.K: array of dimension nxsxd. 
%                     dAKC.C: array of dimension sxnxd. 
%
% REMARK: for numerical reasons in the KPSW data set the transposed system
% is parameterized. 
%
% AUTHOR: dbauer 27.10.2019.

if nargout>3
    np = length(parbull);
    dA = zeros(n,n,np);
    dK = zeros(n,s,np);
    dC = zeros(s,n,np);
end

% if s>=n -> parameters for all matrices contained. 
if s>=n
    % parameters for C: (s-n)*n
    parc= parbull(1:(s-n)*n);
    Cbull = [eye(n);reshape(parc,s-n,n)];
    parbull(1:(s-n)*n)=[];
    % parameters for A:
    Abull = reshape(parbull(1:n^2),n,n);
    parbull(1:n^2)=[];
    
    % parameters for Kbull
    Kbull = reshape(parbull,n,s);
    if nargout>3
        dcp = eye((s-n)*n);
        for j=1:(s-n)*n
            dC(n+[1:(s-n)],:,j)=reshape(dcp(:,j),s-n,n);
        end;
        dap = eye(n^2);
        for j=1:n^2
            dA(:,:,j+(s-n)*n)=reshape(dap(:,j),n,n);
        end;
        dkp = eye(n*s);
        for j=1:n*s
            dK(:,:,s*n+j)=reshape(dkp(:,j),n,s);
        end                   
    end
else
    % if s<n -> no parameter is C. 
    Cbull = [eye(s),zeros(s,n-s)];
    
    % A contains s*n parameters. 
    Abull = zeros(n,n);
    Abull(1:(n-s),s+1:n)=eye(n-s);
    para = parbull(1:n*s);
    Abull(n-s+1:end,:)=reshape(para,s,n);
    
    parbull(1:n*s)=[];
    
    % K contains the remaining parameters. 
    Kbull = reshape(parbull,n,s);
    if nargout>3
        dap = eye(n*s);
        for j=1:n*s
            dA(n-s+1:end,:,j)=reshape(dap(:,j),s,n);
        end;
        dkp = eye(n*s);
        for j=1:n*s
            dK(:,:,s*n+j)=reshape(dkp(:,j),n,s);
        end                   
    end
end;

% switch roles
Abull = Abull';
Cb = Cbull;
Cbull = Kbull';
Kbull = Cb'; 

if nargout>3
    dAKC.A = permute(dA,[2,1,3]);
    dAKC.K = permute(dC,[2,1,3]);
    dAKC.C = permute(dK,[2,1,3]); 
end

function [par,A,K,C] = th2param_I2(th,ds,norm);
% converts th structure into parameter vector
% assumes echelon form for stationary system, heading cxc system to
% correspond to integrated part.
%
% SYNTAX:  par = th2param(th,c);
%
% INPUT:   th ... theta structure
%          ds  ... 2x1 integer vector number of [I(2),I(1)] common trends
%          norm ... indicator: shall the system first be normalized to
%          canform? 
%
% OUTPUT:  param ... kd x 1 real vector containing the parameters.
%
% AUTHOR: dbauer, 7.2.2020.

A = th.A;
K = th.K;
C = th.C; 

[n,s] = size(K);
c = ds(1)*2+ds(2);

nbull = n-c; 

% --- check if normalization necesary? ---
if (nargin < 3) | (norm==1)
    % transform A to Jordan normal form 
    [Tv,A] = Jordan_Au(th.A,ds);
    K = Tv*th.K;
    C = th.C*inv(Tv);
    
    % impose other restrictions in canonical form 
    % start with C1.
    C1 = C(:,1:ds(1));
    [q,r]=qr(C1);
    C(:,1:ds(1)) = q(:,1:ds(1));
    K(1:ds(1),:) = r(1:ds(1),1:ds(1))*K(1:ds(1),:);    
    K(ds(1)+[1:ds(1)],:) = r(1:ds(1),1:ds(1))*K(ds(1)+[1:ds(1)],:);  
    C(:,ds(1)+[1:ds(1)]) = C(:,ds(1)+[1:ds(1)])*inv(r(1:ds(1),1:ds(1))); 
    % next is C3. 
    Tv = eye(n);
    Tv(1:ds(1),2*ds(1)+[1:ds(2)])=-C(:,1:ds(1))'*C(:,2*ds(1)+[1:ds(2)]);
    iTv = inv(Tv);
    A = iTv*A*Tv;
    K = iTv*K;
    C = C*Tv;
    C3 = C(:,2*ds(1)+[1:ds(2)]);
    [q,r]=qr(C3);
    C(:,2*ds(1)+[1:ds(2)]) = q(:,1:ds(2));
    K(2*ds(1)+[1:ds(2)],:) = r(1:ds(2),1:ds(2))*K(2*ds(1)+[1:ds(2)],:);    
    
    % finally C2. 
    Tv = eye(n);
    Tv(2*ds(1)+[1:ds(2)],ds(1)+[1:ds(1)])=-C(:,2*ds(1)+[1:ds(2)])'*C(:,ds(1)+[1:ds(1)]);
    Tv([1:ds(1)],ds(1)+[1:ds(1)])=-C(:,[1:ds(1)])'*C(:,ds(1)+[1:ds(1)]);
    iTv = inv(Tv);
    A = iTv*A*Tv;
    K = iTv*K;
    C = C*Tv;
    
    % now standardize integrated part
    K2 = K(ds(1)+[1:ds(1)],:);
    [Q2,R] = qr(K2);

    K3 = K(2*ds(1)+[1:ds(2)],:);
    [Q3,R] = qr(K3);
    
    Tv =eye(n);
    Tv([1:ds(1)],[1:ds(1)])=Q2;
    Tv(ds(1)+[1:ds(1)],ds(1)+[1:ds(1)])=Q2;
    Tv(2*ds(1)+[1:ds(2)],2*ds(1)+[1:ds(2)])=Q3;
    iTv = inv(Tv);
    A = iTv*A*Tv;
    K = iTv*K;
    C = C*Tv;
    
    % standardize stable part
    Abull = A(c+1:end,c+1:end);
    Kbull = K(c+1:end,:);
    Cbull = C(:,c+1:end);
    
    [~,Ab,Cb,Kb] = ss2ech(Abull',Cbull',Kbull');

    % fill in 
    
    A(c+1:end,c+1:end)=Ab';
    K(c+1:end,:) = Kb';
    C(:,c+1:end) = Cb';
end

% --- convert to real, if necessary ----
A= real(A);
C = real(C);
K = real(K);

% --- now start the parameterization ---
% --- C1 ---
C1 = C(:,1:ds(1));
[parc,QL] = ortho2par_LR(C1);

% --- C3 ---
C3 = QL*C(:,2*ds(1)+[1:ds(2)]); 
% multiplication with QL leads to new matrix starting with zero rows.
[parc3,QL3] = ortho2par_LR(C3(ds(1)+1:end,:));
parc = [parc(:)',parc3(:)'];

% --- C2 ---
C2 = QL*C(:,ds(1)+[1:ds(1)]);
C2h = C2(ds(1)+1:end,:);
C2hh = QL3(ds(2)+1:end,:)*C2h; 
parc = [parc,C2hh(:)'];

% --- K1 ---
K1 = K(1:ds(1),:);

park = K1(:)';

% --- K2 ---
K2 = K(ds(1)+[1:ds(1)],:);
for j=1:size(K2,1)
    park = [park,K2(j,j:end)];
end


% --- K3 ---
K3 = K(2*ds(1)+[1:ds(2)],:);
for j=1:size(K3,1)
    park = [park,K3(j,j:end)];
end


% --- parbull ---
Cbull = C(:,c+1:end);
Kbull = K(c+1:end,:);
Abull = A(c+1:end,c+1:end);
parambull = mat2param_bull(Abull',Cbull',Kbull');

par = [parc(:);park(:);parambull(:)];
par = real(par);


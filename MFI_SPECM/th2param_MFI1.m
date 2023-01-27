function par = th2param_MFI1(th,urs,norm,restrict)
% converts th structure into parameter vector
% assumes echelon form for stationary system, heading cxc system to
% correspond to integrated part.
%
% SYNTAX:  par = th2param(th,c);
%
% INPUT:   th ... theta structure
%          urs  ... unit root structure.
%          norm ... indicator: shall the system first be normalized to
%          canform? 
%          restrict ... contains the specification of restrictions. 
%
% OUTPUT:  param ... kd x 1 real vector containing the parameters.
%         
% REMARK: added changes due to equal Ck matrices. 
%
% AUTHOR: dbauer, 27.10.2019.

if nargin<4
    restrict.det_res = 0;
end

if norm>0
    th = trans_can_form_MFI(th,urs);
end

A = th.A;
K = th.K;
C = th.C; 

[~,s] = size(K);

c = sum((urs(:,2)+1).*urs(:,3));

% --- are we in the restriction of all Cks being equal? 
equalC =0;
if isfield(restrict,'equal')
    if restrict.equal == 1
        equalC = 1;
        if urs(1,2)==0 % real unit root 
            ck = urs(1,3);%ck must be the same then for all k.
            Ckreal = C(:,1:ck);
        else % complex unit root first
            Ckcomp = C(:,1:ck)+ sqrt(-1)*C(:,ck+[1:ck]);
            tCk = real(Ckcomp*inv(Ckcomp(1:ck,:)));
            [q,r]= qr(tCk);
            Ckreal = q(:,1:ck);
        end
        parc = ortho2par_LR(Ckreal);
    end
end

% --- start extracting parameters ----
cur =0;
param1 = [];
if equalC == 1
    param1 = parc(:);
end;

for j=1:size(urs,1) % cycle over unit roots
    ck = urs(j,3); 
    if urs(j,2)>0 %complex root
    % get submatrices 
        Ck = C(:,cur+[1:ck])+ sqrt(-1)*C(:,cur+ck+[1:ck]);
        Kk = (K(cur+[1:ck],:)+ sqrt(-1)*K(cur+ck+[1:ck],:))/2;
    
        % extract parameters from Ck
        if equalC == 0
            parc = ortho2par_LR_c(Ck);
            param1 = [param1;parc(:)];
        else 
            Ckr = Ckreal\Ck;
            parc2 = ortho2par_LR_c(Ckr);
            param1 = [param1;parc2(:)];
        end
           
        % extract parameters from Kk. 
        param1 = [param1;diag(Kk(:,1:ck))];
        for a=1:ck
            for b=a+1:ck
                param1 = [param1;real(Kk(a,b));imag(Kk(a,b))];
            end
        end
        Kk2 = Kk(:,(ck+1):s);
        param1 = [param1(:);real(Kk2(:));imag(Kk2(:))];
        %for a=1:ck
        %    for b=(ck+1):s
        %        param1 = [param1;real(Kk(a,b));imag(Kk(a,b))];
        %    end
        %end
        cur = cur+2*ck;
    else % real root
        Ck = C(:,cur+[1:ck]);
        Kk = K(cur+[1:ck],:);
    
        % extract parameters from Ck
        if equalC == 0
            parc = ortho2par_LR(Ck);
            param1 = [param1;parc(:)];
        else
            Ckq = Ckreal\Ck;
            Kk = Ckq*Kk;
        end
        % extract parameters from Kk. 
        param1 = [param1;diag(Kk(:,1:ck))];
        for a=1:ck
            for b=a+1:ck
                param1 = [param1;real(Kk(a,b))];
            end
        end
        for a=1:ck
            for b=(ck+1):s
                param1 = [param1;real(Kk(a,b))];
            end
        end
        cur =cur+ck;
    end
end


% --- obtain stable subsystem ---
Cbull = C(:,c+1:end);
Kbull = K(c+1:end,:);
Abull = A(c+1:end,c+1:end);

% --- convert to parameters ---
parambull = mat2param_bull(Abull',Cbull',Kbull');

paromi = extr_lowtri(th.Omega);

par = [paromi(:);param1(:);parambull(:)];
par = real(par);

% add parameters for D. 
par = [par;th.D(:)];



function th = trans_can_form_MFI(th,urs);
% trans_can_form_MFI transforms the system th to the canonical form of Bauer and Wagner (2012).
% 
% SYNTAX: th = trans_can_form_MFI(th,urs);
%
% INPUT: th ... theta_urs structure.
%        urs ...  M x 3 matrix containing the unit root structure.
%        [frequency, complex valued?, number of common trends].
%
% OUTPUT: th ... transformed state space system.
%
% AUTHOR: dbauer, 9.10.2020. 

if th.which == "poly"
    th = poly2ss(th);
end;

A = th.A;
K = th.K;
C = th.C; 
n = size(A,1);

[v,d]= seig(A);
% A = v*d*inv(v)
% eigenvalues ordered decreasing in absolute value. 
se = diag(d);
ase = abs(se);



in = find(ase>0.999); 
stab = setdiff(1:n,in);

% check, if enough. 
if length(in)~= sum((urs(:,2)+1).*urs(:,3))
    warning('Unit root structure and location of roots does not conform -> ignoring unit root structure.');
end;

% first separate unit roots from rest. 
T = v(:,[in(:);stab(:)]);
A = inv(T)*A*(T);
K = inv(T)*K;
C = C*(T);

c = length(in); 
% order them in the right order
se = diag(A);


in2= [];
for j=1:size(urs,1)
    zk = exp(sqrt(-1)*urs(j,1)*pi);
    if urs(j,2)>0 % complex root
        ind1 = find(abs(se-zk)<10^(-4));
        ind2 = find(abs(se-conj(zk))<10^(-4));
        in2 = [in2(:);ind1(:);ind2(:)];
    else
      ind = find(abs(se-zk)<10^(-4));
      in2 = [in2(:);ind(:)];
    end
end

stab = setdiff(1:n,in2);
bind = [in2;stab(:)];
c = length(in2);
A = A(bind,bind);
K = K(bind,:);
C = C(:,bind); 

% now iterate over unit roots.
se = diag(A);
cur = 0;

for j=1:size(urs,1)
    zk = exp(sqrt(-1)*urs(j,1)*pi);
    if urs(j,2)>0 % complex root
        ind1 = find(abs(se-zk)<10^(-4));
        ind2 = find(abs(se-conj(zk))<10^(-4));
        ind = [ind1,ind2]; 
        out = setdiff(1:n,[ind1,ind2]);
        %correct A
        A(ind,ind)=0;
        A(ind1,ind1)= real(zk)*eye(length(ind1));
        A(ind2,ind2)= real(zk)*eye(length(ind1));
        A(ind1,ind2)= imag(zk)*eye(length(ind1));
        A(ind2,ind1)= -imag(zk)*eye(length(ind1));
        
        % transform C: 
        [q,r]=qr(C(:,ind1));
        C(:,ind1)= q(:,1:length(ind1)); 
        K(ind1,:)=r(1:length(ind1),1:length(ind1))*K(ind1,:);
        % correct K(ind,:);
        [q,r]=qr(K(ind1,1:length(ind1)));
        C(:,ind1)= C(:,ind1)*q; 
        K(ind1,:)=inv(q)*K(ind1,:);
        
        % fill in complex conjugate part
        C(:,ind2)=imag(C(:,ind1));
        C(:,ind1)=real(C(:,ind1));

        K(ind2,:) = -2*imag(K(ind1,:)); 
        K(ind1,:) = 2*real(K(ind1,:)); 
                
        % finally move to the front.
        %bind = [ind(:);out(:)];
        %A = A(bind,bind);
        %K = K(bind,:);
        %C = C(:,bind);
        %se = se(bind);
    else % real unit root 
        ind = find(abs(se-zk)<10^(-4));
        out = find(abs(se-zk)>10^(-4));
        %correct A
        A(ind,ind)= real(zk)*eye(length(ind));
        % transform C: 
        [q,r]=qr(C(:,ind));
        C(:,ind)= q(:,1:length(ind)); 
        K(ind,:)=r(1:length(ind),1:length(ind))*K(ind,:);
        % correct K(ind,:);
        [q,r]=qr(K(ind,1:length(ind)));
        C(:,ind)= real(C(:,ind)*q); 
        K(ind,:)=real(inv(q)*K(ind,:));
        
        % finally move to the front.  
        %A = A([ind;out],[ind;out]);
        %K = K([ind;out],:);
        %C = C(:,[ind;out]);
        %se = se([ind;out]);
    end
end


% and finally the remaining stable part. 
stab = [c+1:n];
Abull = A(stab,stab);

rhoa = max(abs(eig(Abull)));
if rhoa>1.00001
    warning("System is explosive. Adjusting!");
    Abull = Abull/rhoa*0.99;
end;



Kbull = K(stab,:);
Cbull = C(:,stab); 

[~,Ab,Kb,Cb] = ss2ech_n(Abull',Cbull',Kbull'); 

A(stab,stab) = Ab';
K(stab,:)= Cb';
C(:,stab) = Kb'; 

th.A = real(A);
th.K = real(K);
th.C = real(C); 



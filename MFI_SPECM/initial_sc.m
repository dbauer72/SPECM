function estimate = initial_sc(pX,S,k,pXk,Xm2,urs,Joh_j);
% provides the initial estimate of the seasonal coint vecm representation
% of Johansen and Schaumburg.
%
% SYNTAX: estimate = initial_sc(pX,pXk,Xm2,urs);
%
% INPUT: pX  ... T x s  ... p(L)X_t
%        S  ... integer; number of observations per year
%        k  ... integer; number of lagged differences
%        pXk ... T x ks ... p(L)X_{t-1}, ..., p(L)X_{t-k}
%        Xm2 ... T x k x s ... X_t^{(m)}
%        urs ... unit root structure.
%        Joh_j ... specification of deterministics.
%
% OUTPUT: estimate ... vecm_mfi structure.
%
% AUTHOR: dbauer, 27.11.2015.

% initialize structure
estimate = vecm_mfi();
estimate.urs = urs; 
estimate.s = size(pX,2);
estimate.k = k;

% full set of frequencies 
om = [0:(S/2)]'/S*2; % frequencies between [0,1].
fr = exp(om*pi*sqrt(-1)); % corr. unit roots. 

ursf(:,1)= om(:);
ursf(:,2)= 1;
ursf([1,end],2)=0;

% obtain unrestricted estimate
[T,zs,MS]= size(Xm2);
s = estimate.s;
ind = [size(pXk,2)];
X = pXk;

M = size(urs,1);
ursf(:,3)=0;
for m=1:M
    ursf(urs(m,4),3)=m;
end;
ind = cell(MS,1);

cur_omit = k*s; % number of differenced estimators. 

for m=1:MS % cycle over roots
    if ursf(m,3)>0 % frequency contained in urs
        xx = squeeze(Xm2(:,:,m));
        curind =size(X,2)+1;
        if ursf(m,2)>0 % complex root 
            X = [X,xx];        
        else 
            X = [X,xx(:,1:s)];
            if Joh_j == 4
                X = [X,xx(:,2*s+1)];
            end
            if (Joh_j < 4)&&(ursf(m,1)>0.00001)
                X = [X,xx(:,2*s+1)];
            end
        end;
        ind{m} = [curind:size(X,2)]; 
    else % frequency not contained -> estimate without restriction. 
        if ursf(m,2)>0
            if Joh_j<5
                plus1 = 2;
            else
                plus1 = 0;
            end;
            ind{m} = cur_omit+[1:2*s+plus1];
            cur_omit = cur_omit+2*s+plus1;
        else % real root
            if Joh_j <5
                plus1 = 1;
            else
                plus1 = 0;
            end;
            
            if (Joh_j < 4)&&(m==1) 
                plus1=0;
            end;
            ind{m} = cur_omit+[1:s+plus1];
            cur_omit = cur_omit+s+plus1;
        end
    end
end;

bhat = (X\pX)';
estimate.Gamma = bhat(:,1:k*s);


for m=1:MS  % cycle over roots
    if ursf(m,3)>0
        
        cm = s-urs(ursf(m,3),3);
        if ursf(m,2)>0 % complex root 
            % set up complex valued problem 
            out = setdiff(1:size(X,2),ind{m}); 
            R0t = pX - X(:,out)*(X(:,out)\pX);
            R1t = X(:,ind{m}) - X(:,out)*(X(:,out)\X(:,ind{m}));
            Rit = R1t(:,1:s)+sqrt(-1)*R1t(:,s+[1:s]); % complex valued X_t^{(m)} 
            if size(R1t,2)>2*s % determinstics contained?
                Rit = [Rit,R1t(:,2*s+1)+sqrt(-1)*R1t(:,2*s+2)];
            end;
            S01 = R0t'*Rit/(T-k);
            S11 = Rit'*Rit/(T-k);
            bb = S01*inv(S11);
            % reducde rank of this estimator. 
            [U,S,V] = svd(bb);
            ah =  U(:,1:cm)*S(1:cm,1:cm);
            bh = V(:,1:cm)';
            bb = bh(:,1:cm);
            bh = inv(bb)*bh;
            ah = ah*bb;
            % now fill in result. 
            br = real(bh');
            bi = imag(bh');
            bhat(:,ind{m})= 2*[real(ah),imag(ah)]*[br,-bi;bi,br]';
            estimate.alpha{m} =2*[real(ah),imag(ah)];
            estimate.beta{m} = [br,-bi;bi,br];
        else
            alphah = bhat(:,ind{m});
            [U,S,V] = svd(alphah);
            ah =  U(:,1:cm)*S(1:cm,1:cm);
            bh = V(:,1:cm)';
            bb = bh(:,1:cm);
            bh = inv(bb)*bh;
            ah = ah*bb;
            estimate.alpha{m} =ah;
            estimate.beta{m} = bh';
            bhat(:,ind{m})= ah*bh;            
        end;
    else % frequency not contained -> keep unrestricted estimate.
        estimate.alpha{m} =bhat(:,ind{m});
        estimate.beta{m} = eye(length(ind{m}));
    end;
end

u = pX - X * bhat';
estimate.Omega = u'*u/T;


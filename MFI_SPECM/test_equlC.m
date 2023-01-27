% 
M = 100;
T = 2000;
% test estimation of system with equal Ck.
n=5;
s=2;
A = diag([1,-1,sqrt(-1),-sqrt(-1),0.5*eye(n-4)]);
K = randn(n,2);
Kc = K(3,:)+sqrt(-1)*K(4,:);
K(3,:)=(Kc);
K(4,:)=conj(Kc);
K(2,1)= abs(K(2,1)); % turn 2,1 element positive to conform with canonical form. 

om = randn(1,2);
om = om/norm(om);
C = [[1;1]*[1,1,om],randn(2,n-4)];
Cc = C(:,3)+sqrt(-1)*C(:,4);
C(:,4)=(Cc);
C(:,3)= conj(Cc); 

Tr = [eye(2),zeros(2,n-2);[zeros(2,2),[1,1;-sqrt(-1),sqrt(-1)],zeros(2,n-4)];zeros(n-4,4),eye(n-4)];
Ar = Tr*A*inv(Tr);
Kr = Tr*K
Cr = C*inv(Tr);

th = theta_urs();
th.A = Ar;
th.K = Kr;
th.C = Cr;
th.Omega= eye(s); 

urs = [0,0,1;0.5,1,1;1,0,1];
th.urs = urs; 
par = th2param_MFI1(th,urs,0);

thn = param2th_MFI1(par,s,n,urs);

restrict.det_res= 0;

% loop over replications 

load invertible_theta 

h = waitbar(0,'Please wait ...');
for m=1:M
    waitbar(m/M,h);
    % generate y
    x = zeros(n,1);
    A= thn.A
    K=thn.K;
    C = thn.C;

    
    y = zeros(T,s);
    for t=1:T
        e=randn(s,1);
        y(t,:)= x'*C'+e';
        x = A*x + K*e;
    end

    % --- estimate without restriction ---
    [result_ur,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM_MFI1(y,s,4,n,urs,4,0,restrict);

    if (result_ur.deviance> 4000) % most likely cause: no invertible system can be found. 
        [result_ur,thc,Ac,Kc,Cc,Omegac,thi,thi2,lle] = SPECM_MFI1(y,s,4,n,urs,4,0,restrict, th_repl);
    end

    % --- estimate with restriction ---
    restrict.equal = 1; 
    [result_r,thr,Ar,Kr,Cr,Omegar,thi_approx] = SPECM_MFI1_equalC(y,s,4,n,urs,4,0,restrict,result_ur.theta);
    % --- reestimate from here, if not successfull ---
    restrict.equal = 0;
    [result_ur2,thc2,Ac2,Kc2,Cc2,Omegac2,lle2] = SPECM_MFI1_freeC(y,s,4,n,urs,4,0,restrict,thr);
    
    it = 0;
    while ((abs(result_ur.deviance-result_ur2.deviance)>0.001)&&(it<15)) 
        result_ur = result_ur2;
        restrict.equal = 1; 
        [result_r,thr,Ar,Kr,Cr,Omegar,thi_approx] = SPECM_MFI1_equalC(y,s,4,n,urs,4,0,restrict,result_ur.theta);
        % --- reestimate from here, if not successfull ---
        restrict.equal = 0;
        [result_ur2,thc2,Ac2,Kc2,Cc2,Omegac2,lle2] = SPECM_MFI1_freeC(y,s,4,n,urs,4,0,restrict,thr);
        it = it+1;
    end

    % compare results and use better one. 
    dev_ur = min(result_ur.deviance,result_ur2.deviance);
    dev_r = result_r.deviance;
    it2 = 0;
    while (dev_ur>dev_r)&&(it2<5)
        result_ur2 = random_improve_call(result_ur2,0.0001,15);
        dev_ur = min(dev_ur,result_ur2.deviance);
        dev_r = result_r.deviance;    
        it2 = it2+1;
    end

    devs(m,:)= [dev_r,dev_ur,it,it2];

end
close(h);

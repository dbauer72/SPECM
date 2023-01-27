
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test testing in hypotheses on beta  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T= 5000;
p=4;
A = [1,0,0,0;0,1,0,0;0,0,.5,0;0,0,0,.5];
C(1:p,1) = [1,0,0,0]';
C(:,2) = [0,1,0,0]';
C(:,3:4) = randn(4,2);
K(1,1:4) = [1,0,0,0];
K(2,1:4) = [0,1,0,0];
K(3:4,1:4) = randn(2,4);

while max(abs(eig(A-K*C)))>0.95
    K = K*.5;
end

Abar = A-K*C;
r=2;
rest = 'y';
Joh_j = 5;
restrict.beta = [0,0;0,0;1,0;0,1];
restrict.type = 'A';

% test drive
clear u y;
u = randn(T,p);
x = randn(1,p);
y(1,1:p)=u(1,:);
for t=2:T
    x = x*A' + u(t-1,:)*K';
    y(t,:) = x*C'+ u(t,:);
end

[LL,alphahat,betahat] = RH_specm_H0(y,Abar,K,r,rest,Joh_j,restrict);

% now start the engine. 
p=4;
M = 1000;
T=5000;
Joh_j =5;

restrict.beta = [0,0;0,0;1,0;0,1];
restrict.type = 'A';


for m = 1:M
    u = randn(T,p);
    x = randn(1,p);
    y(1,1:p)=u(1,:);
    for t=2:T
        x = x*A' + u(t-1,:)*K';
        y(t,:) = x*C'+ u(t,:);
    end

    [LL(m,:),~,~,df] = RH_specm_H0(y,Abar,K,r,rest,Joh_j,restrict);
end

% distribution of difft(LL)?
dL = LL(:,2)-LL(:,1);
kdfft1(dL,'knorm',1024,.2,1); hold on;

df = r*(p-r);
fr = 0:0.1:20;
chi2 = chis_pdf(fr,df);
plot(fr,chi2);
set(gca,'xlim',[0,max(dL)]);

%%% changes for Joh_j =2, since then beta needs to contain another row. %%
p=4;
M = 1000;
T=5000;
Joh_j =4;
restrict.beta = [0,0;0,0;0,0;1,0;0,1];
restrict.type = 'A';

for m = 1:M
    u = randn(T,p);
    x = randn(1,p);
    y(1,1:p)=u(1,:);
    for t=2:T
        x = x*A' + u(t-1,:)*K';
        y(t,:) = x*C'+ u(t,:);
    end

    [LL(m,:),~,~,df] = RH_specm_H0(y,Abar,K,r,rest,Joh_j,restrict);
end

% distribution of difft(LL)?
dL = LL(:,2)-LL(:,1);
kdfft1(dL,'knorm',1024,.1,1); hold on;

%df = r*(p-r+1);
fr = 0:0.1:20;
chi2 = chis_pdf(fr,df);
plot(fr,chi2);
set(gca,'xlim',[0,max(dL)]);

%%%% type B.
restrict.type = 'B';
Joh_j = 5;
restrict.H = [zeros(1,3);eye(3)];

for m = 1:M
    u = randn(T,p);
    x = randn(1,p);
    y(1,1:p)=u(1,:);
    for t=2:T
        x = x*A' + u(t-1,:)*K';
        y(t,:) = x*C'+ u(t,:);
    end

    [LL(m,:),~,~,df] = RH_specm_H0(y,Abar,K,r,rest,Joh_j,restrict);
end

% distribution of difft(LL)?
dL = LL(:,2)-LL(:,1);
kdfft1(dL,'knorm',1024,.1,1); hold on;
s = size(restrict.H,2);
%df = r*(p-s);
fr = 0:0.1:20;
chi2 = chis_pdf(fr,df);
plot(fr,chi2);
set(gca,'xlim',[0,max(dL)]);


%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%% type C   %%%%
%%%%%%%%%%%%%%%%%%

restrict.type = 'C';
Joh_j = 2;
restrict.b = [0,0,0,1,0]';

for m = 1:M
    u = randn(T,p);
    x = randn(1,p);
    y(1,1:p)=u(1,:);
    for t=2:T
        x = x*A' + u(t-1,:)*K';
        y(t,:) = x*C'+ u(t,:);
    end

    [LL(m,:),~,~,df] = RH_specm_H0(y,Abar,K,r,rest,Joh_j,restrict);
end

% distribution of difft(LL)?
dL = LL(:,2)-LL(:,1);
kdfft1(dL,'knorm',1024,.1,1); hold on;
s = size(restrict.b,2);
%df = s*(p-r+1);
fr = 0:0.1:20;
chi2 = chis_pdf(fr,df);
plot(fr,chi2);
set(gca,'xlim',[0,max(dL)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%% tests for alpha          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=eye(p);
Pi = I-C*inv(I-Abar)*K

alpha_0 = [zeros(2,2);eye(2)];

restrict.A = [zeros(1,3);eye(3)];
restrict.type = 'A'
clear LL
Joh_j=5;

for m = 1:M
    u = randn(T,p);
    x = randn(1,p);
    y(1,1:p)=u(1,:);
    for t=2:T
        x = x*A' + u(t-1,:)*K';
        y(t,:) = x*C'+ u(t,:);
    end

    [LL(m,:),~,~,df] = RH_specm_Ha(y,Abar,K,r,rest,Joh_j,restrict);
end

% distribution of difft(LL)?
dL = LL(:,2)-LL(:,1);
kdfft1(dL,'knorm',1024,.1,1); hold on;
fr = 0:0.1:20;
chi2 = chis_pdf(fr,df);
plot(fr,chi2);
set(gca,'xlim',[0,max(dL)]);

%%%%% type B %%%%%%
restrict.type = 'B';
restrict.a = [0,0,1,0]';

Joh_j=4;

for m = 1:M
    u = randn(T,p);
    x = randn(1,p);
    y(1,1:p)=u(1,:);
    for t=2:T
        x = x*A' + u(t-1,:)*K';
        y(t,:) = x*C'+ u(t,:);
    end

    [LL(m,:),~,~,df] = RH_specm_Ha(y,Abar,K,r,rest,Joh_j,restrict);
end

% distribution of difft(LL)?
dL = LL(:,2)-LL(:,1);
kdfft1(dL,'knorm',1024,.1,1); hold on;
s = size(restrict.b,2);
fr = 0:0.1:20;
chi2 = chis_pdf(fr,df);
plot(fr,chi2);
set(gca,'xlim',[0,max(dL)]);


%%%%% type C %%%%%%
I=eye(p);
Pi = I-C*inv(I-Abar)*K

restrict.type = 'C'
restrict.A = [0,0;0,0;1,0;0,1];
Joh_j=5;

T= 5000;
M=10000;
for m = 1:M
    u = randn(T,p);
    x = randn(1,p);
    y(1,1:p)=u(1,:);
    for t=2:T
        x = x*A' + u(t-1,:)*K';
        y(t,:) = x*C'+ u(t,:);
    end

    [LL(m,:),~,~,df] = RH_specm_Ha(y,Abar,K,r,rest,Joh_j,restrict);
end

% distribution of difft(LL)?
dL = LL(:,2)-LL(:,1);
kdfft1(dL,'knorm',1024,.1,1); hold on;
fr = 0:0.1:20;
chi2 = chis_pdf(fr,df);
plot(fr,chi2);
set(gca,'xlim',[0,max(dL)]);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% testing hypotheses of real beta or alpha   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restrict.type = 'A';
restrict.betafr = 0.5;
restrict.alphafr = 0.5; 


A = [1,0,0,0;0,-1,0,0;0,0,0,1;0,0,-1,0];
K= randn(4,4);
C = eye(4);
C(:,4) =[0,0,-.1,0]'; 

while max(abs(eig(A-K*C)))>0.99
    K = randn(4,4)*.5;
end;

T= 2000;
p = 4;

u = randn(T,p);
x = randn(1,p);
y(1,1:p)=u(1,:);
for t=2:T
    x = x*A' + u(t-1,:)*K';
    y(t,:) = x*C'+ u(t,:);
end

urs = [0,0,1;.5,1,1;1,0,1];

restrict.type = 'A';
[LLA] = RH_SPECM_MFI1_H(y,4,A-K*C,K,Joh_j,urs,restrict);

restrict.type = 'B';
[LLB] = RH_SPECM_MFI1_H(y,4,A-K*C,K,Joh_j,urs,restrict);

restrict.type = 'C';
[LLC] = RH_SPECM_MFI1_H(y,4,A-K*C,K,Joh_j,urs,restrict);

% now test on more samples 
M=10000;
LL = zeros(M,4);
T= 5000;

% type 'A'
for m=1:M
    u = randn(T,p);
    x = randn(1,p);
    y(1,1:p)=u(1,:);
    for t=2:T
        x = x*A' + u(t-1,:)*K';
        y(t,:) = x*C'+ u(t,:);
    end

    restrict.type = 'A';
    [LLA,df] = RH_SPECM_MFI1_H(y,4,A-K*C,K,Joh_j,urs,restrict);

    [LL0] = RH_SPECM_MFI1(y,4,A-K*C,K,Joh_j);
    LL(m,:)= [LL0{2}(4),LLA{2},LLB{2},LLC{2}];
end

df = df(2);
dL = -LL(:,2)+LL(:,1);
kdfft1(dL,'knorm',1024,.1,1); hold on;
fr = 0:0.1:20;
chi2 = chis_pdf(fr,df);
plot(fr,chi2);
set(gca,'xlim',[0,max(dL)]);

% type 'B' 
A = [1,0,0,0;0,-1,0,0;0,0,0,1;0,0,-1,0];
K= randn(4,4);
C = eye(4);
C(:,4) =[0,0,0,0]'; 

while max(abs(eig(A-K*C)))>0.99
    K = randn(4,4)*.5;
    K(4,:)=0;
end;

Abar = A-K*C;
z=sqrt(-1);
Pi = I - C*inv(I-z*Abar)*K*z

LL = zeros(M,2);

% type 'A'
for m=1:M
    u = randn(T,p);
    x = randn(1,p);
    y(1,1:p)=u(1,:);
    for t=2:T
        x = x*A' + u(t-1,:)*K';
        y(t,:) = x*C'+ u(t,:);
    end

    restrict.type = 'B';
    [LLA,df] = RH_SPECM_MFI1_H(y,4,A-K*C,K,Joh_j,urs,restrict);

    [LL0] = RH_SPECM_MFI1(y,4,A-K*C,K,Joh_j);
    LL(m,:)= [LL0{2}(4),LLA{2}];
end

df = df(2);
dL = -LL(:,2)+LL(:,1);
kdfft1(dL,'knorm',1024,.1,1); hold on;
fr = 0:0.1:20;
chi2 = chis_pdf(fr,df);
plot(fr,chi2);
set(gca,'xlim',[0,max(dL)]);

%%% type 'C' %%%
LL = zeros(M,2);

% type 'A'
for m=1:M
    u = randn(T,p);
    x = randn(1,p);
    y(1,1:p)=u(1,:);
    for t=2:T
        x = x*A' + u(t-1,:)*K';
        y(t,:) = x*C'+ u(t,:);
    end

    restrict.type = 'C';
    [LLA,df] = RH_SPECM_MFI1_H(y,4,A-K*C,K,Joh_j,urs,restrict);

    [LL0] = RH_SPECM_MFI1(y,4,A-K*C,K,Joh_j);
    LL(m,:)= [LL0{2}(4),LLA{2}];
end

df = df(2);
dL = -LL(:,2)+LL(:,1);
kdfft1(dL,'knorm',1024,.1,1); hold on;
fr = 0:0.1:20;
chi2 = chis_pdf(fr,df);
plot(fr,chi2);
set(gca,'xlim',[0,max(dL)]);

%%% test imposition of restrictions at z=1. 
restrict.type = 'I1.1.A';
restrict.beta = [zeros(1,3);eye(3)];
[LLA,df] = RH_SPECM_MFI1_H(y,4,A-K*C,K,Joh_j,urs,restrict);

LL = zeros(M,2);

% type 'I1.1.A'
for m=1:M
    u = randn(T,p);
    x = randn(1,p);
    y(1,1:p)=u(1,:);
    for t=2:T
        x = x*A' + u(t-1,:)*K';
        y(t,:) = x*C'+ u(t,:);
    end

    [LLA,df] = RH_SPECM_MFI1_H(y,4,A-K*C,K,Joh_j,urs,restrict);

    [LL0] = RH_SPECM_MFI1(y,4,A-K*C,K,Joh_j);
    LL(m,:)= [LL0{2}(4),diff(LLA{1})];
end

df = df(1);
dL = LL(:,2);
kdfft1(dL,'knorm',1024,.1,1); hold on;
fr = 0:0.1:20;
chi2 = chis_pdf(fr,df);
plot(fr,chi2);
set(gca,'xlim',[0,max(dL)]);

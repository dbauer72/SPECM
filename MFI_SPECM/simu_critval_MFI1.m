function [t_val]= simu_critval_MFI1(T,M,om,n,det_fr);
% simulates the critical values for the Lambda test on the number of unit
% roots. 
%

zj = exp(sqrt(-1)*om*pi*2);

if det_fr
    dt = (conj(zj).^[1:T])';
end

for m=1:M
    u = randn(T,n);
    x = filter(1,[1,-zj],u,0);
    if det_fr
        x = x - dt*(dt\x);
    end;
    t_c(m)=trace(x(1:end-1,:)\u(2:end,:));
end;

t_val =T*abs(t_c);
ecdf(t_val);

%Q = [cos(2*pi*om),sin(2*pi*om);-sin(2*pi*om),cos(2*pi*om)];
%
%for m=1:M
%    u = randn(T,n*2);
%    x = ltitr(kron(Q,eye(n)),eye(2*n),u);
    
%    A =x(1:end-1,:)\x(2:end,:);
%    dlambda = eig(A)-zj;
%    dlambda = dlambda(abs(dlambda)<0.5);
%    t_r(m) = sum(dlambda);
%end;
%
%[Fr,Xr] = ecdf(T*abs(t_r));
%
hold on;
%plot(Xc,Fc);
%plot(Xr,Fr,'r');
function vals = simu_crit(p,nd,T,M,Joh_j);
% simulates the critical values for complex valued trace tests.
%
% SYNTAX: vals = simu_crit(p,dim,T,M,Joh_j);
%
% INPUT:   p  ... P vector of reals; probabilities to estimate distribution
%          nd ... integer; maximal dimension. 
%          T  ... integer; number of sample to approximate Brownian motion
%          M  ... integer; number of replications. 
%          Joh_j ... integer; chooses the deterministics specs.
%
% OUTPUT: vals ... P x nd matrix of quantiles.
%
% AUTHOR: dbauer, 24.8.2020. 

P = length(p);

for s=1:nd
   for m=1:M
       u = randn(T,2*s);
       y = cumsum([zeros(1,2*s);u(1:end-1,:)]);
       switch Joh_j 
           case {1,2,3} 
               my = mean(y); 
               y = y- ones(T,1)*my;
               my = mean(y);
               mu = mean(u);
               uBd = u(:,1:s)'*y(:,1:s) + u(:,s+[1:s])'*y(:,s+[1:s]);
               uBu = -u(:,1:s)'*y(:,s+[1:s]) + u(:,s+[1:s])'*y(:,[1:s]);
               uB = [ uBd, - uBu, mu(1:s)', - mu(s+[1:s])'; uBu, uBd , mu(s+[1:s])' , mu(1:s)'];
               BBd = y(:,1:s)'*y(:,1:s) + y(:,s+[1:s])'*y(:,s+[1:s]);
               BBu = -y(:,1:s)'*y(:,s+[1:s]) + y(:,s+[1:s])'*y(:,[1:s]);
               BB = [ BBd, - BBu , my(1:s)', - my(s+[1:s])' ; BBu,BBd , my(s+[1:s])' ,  my(1:s)' ; my , 1, 0 ; -my(s+[1:s]),my(1:s),0,1];
               vv(m)= 0.5*trace(uB*inv(BB)*uB');
           case 4
               uBd = u(:,1:s)'*y(:,1:s) + u(:,s+[1:s])'*y(:,s+[1:s]);
               uBu = -u(:,1:s)'*y(:,s+[1:s]) + u(:,s+[1:s])'*y(:,[1:s]);
               mu = mean(u);
               uB = [ uBd, - uBu, mu(1:s)', - mu(s+[1:s])'; uBu, uBd , mu(s+[1:s])' , mu(1:s)'];
               BBd = y(:,1:s)'*y(:,1:s) + y(:,s+[1:s])'*y(:,s+[1:s]);
               BBu = -y(:,1:s)'*y(:,s+[1:s]) + y(:,s+[1:s])'*y(:,[1:s]);
               my = mean(y); 
               BB = [ BBd, - BBu , my(1:s)', - my(s+[1:s])' ; BBu,BBd , my(s+[1:s])' ,  my(1:s)' ; my , 1, 0 ; -my(s+[1:s]),my(1:s),0,1];
               vv(m)= 0.5*trace(uB*inv(BB)*uB');
           case 5
               uBd = u(:,1:s)'*y(:,1:s) + u(:,s+[1:s])'*y(:,s+[1:s]);
               uBu = -u(:,1:s)'*y(:,s+[1:s]) + u(:,s+[1:s])'*y(:,[1:s]);
               uB = [ uBd, - uBu; uBu,uBd];
               BBd = y(:,1:s)'*y(:,1:s) + y(:,s+[1:s])'*y(:,s+[1:s]);
               BBu = -y(:,1:s)'*y(:,s+[1:s]) + y(:,s+[1:s])'*y(:,[1:s]);
               BB = [ BBd, - BBu; BBu,BBd];
               vv(m)= 0.5*trace(uB*inv(BB)*uB');
       end
       
   end
    vals(s,:)= prctile(vv,p*100);
end

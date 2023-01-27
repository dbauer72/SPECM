function [LLy,LLn,b_eta,alphahat,betahat,Chat] = test_Johansen_RH_I1(result,alpha,Joh_j,plots);
% calculates the Johansen trace test using the SS-ECM representation of Ribarits
% and Hanzon. 
%
% SYNTAX: [LLy,LLn,b_eta,alphahat,betahat,Chat] = test_Johansen_RH(result,alpha,Joh_j,plots);
%
% INPUT: result ... structure containing the estimation results.
%        alpha  ... real; significance level. Default: alpha = 0.05.
%        Joh_j  ... integer; one of the five Johansen levels for
%                     deterministics.
%        plots  ... indicator if plotting is wanted. 
%
% OUTPUT:  LLy  ... (s+1) x 1 vector of test statistics of H(j) versus H(s).
%          LLn  ... (s+1) x 1 vector of test stats (without restriction in
%                       RH procedure).
%          b_eta ... estimated cointegrating rank
%          alphahat ... corresponding estimate of Pi = alphahat * betahat'
%          betahat  ... cointegrating vectors.
%          Chat     ... estimate of matrix C according to RH procs. 
%
% REMARK: implements the Ribarits-Hanzon version of the Johansen trace
% tests. The two version include the restriction on the matrix C or an
% unrestricted estimate. 
%
% AUTHOR: dbauer, 29.11.2019

thc = result.theta; 
Abar = thc.A-thc.K*thc.C;
B =thc.K;
s = size(thc.C,1); 

tily = result.y; %- result.dt*thc.D'; 

for r=0:s,
    [LLy(r+1),alphahat,betahat,Chat] = RH_specm(tily,Abar,B,r,'y',Joh_j); %,result.restrict);
    [LLn(r+1),alphahat,betahat,Chat] = RH_specm(tily,Abar,B,r,'n',Joh_j); %,result.restrict);
end;

% Table really corresponding to model without constant
load eta_critical % MAT file contains simulated critical values. 
etamat = eta_est{Joh_j}; 
% [248.771 256.231 262.694 270.466 
%          207.212 214.069 220.000 226.945
%		  169.302 175.472 181.199 187.821
%          135.164 140.736 145.797 151.696 
%          104.933 109.932 114.241 119.578 			   
%		  78.303 82.614 86.363 91.115 
%          55.54 59.235 62.68 66.705
%		  36.579 39.71 42.588 45.995 
%          21.581 24.076 26.422 29.194
%		  10.347 12.212 13.939 16.158
%          2.98 4.141 5.302 7.015];
      
% Selection of column corresponding to ALPHA
if alpha == 0.10
    col = 1;
elseif alpha == 0.05
    col = 2;
elseif alpha == 0.025
    col = 3;
elseif alpha == 0.01
    col = 4;
end;

% Testing sequence
% TRACE TEST
T = size(tily,1);
tracevec = (LLy-LLy(end));
tracevecn = (LLn-LLn(end));
i = 1;
while i <= s,
    if i == 1,
        if tracevec(1) <= etamat(11-s+i,col),
            b_eta = 0;
            i = s+1;
        end;
    end;
    if i <= s,
    if tracevec(i) > etamat(11-s+i,col),
        i = i + 1;
    else 
        b_eta = i-1;
        i = s+1;
    end;    
    if i == s,
        if tracevec(i) > etamat(11-s+i,col),
            b_eta = s;
            i = s+1;
        end;            
    end;    
    end;
end;



if plots
    figure;
    hold on;
    set(gca,'XTick',0:s-1);
    bar(0:s-1,etamat(11-s+1:11,col),'b');
    plot(0:s-1,tracevec(1:end-1),'r*');
    title(sprintf('Estimated coint rank: %d',b_eta));
end;

% perform estimation again with c= s-b_eta. 
[~,alphahat,betahat,Chat] = RH_specm(tily,Abar,B,b_eta,'y',Joh_j);



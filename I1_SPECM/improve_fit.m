function result_ur = improve_fit(result,restrict);
% improves the fit starting from new initial point. 
%
% SYNTAX: result_ur = improve_fit(result,restrict);
%
% INPUT: result ... result structure
%        restrict ... provides information on det_res; rest is ignored.
%
% OUTPUT:  result_ur ... updated result structure without restrictions.
%
% REMARK: can be used in order to test the validity of restrictions by
%        taking them away and reestimating. 
%
% AUTHOR: dbauer 14.1.2020. 

rest.det_res = restrict.det_res; 
try 
    [par] = syst2param(result.theta,result.c,rest);
catch 
    disp('improve_fit: Calculations of parameter vector failed. Potentially theta does not match.');
end
% fill in variables 
z = [result.y,result.dt];
[s,n]=size(result.theta.C);
m = size(result.theta.D,2);


% improve estimate
options = optimoptions('fminunc','display','iter');

parc = fminunc(@(x) cal_quasi_like(x,z,s,m,n,result.urs,result.Pbull,rest),par,options);
[llc_r,resc_r]=cal_quasi_like(parc,z,s,m,n,result.urs,result.Pbull,rest);
[A,K,C,D,Omega,thc_r] = param2syst(parc,s,n,m,result.c,rest); 

% finish up 
result_ur = result;
result_ur.theta = thc_r;
result_ur.deviance = llc_r;
result_ur.res = resc_r;
result_ur.restrict = restrict;
result_ur.call = [result.call,';improve_fit;'];
result_ur.param = parc;

function th = convert_ni_theta(nML,c);
% converts the output of Lukas functions to the theta format. 

th = theta();
th.which = 'ss';
th.A = nML.A;
th.K = nML.K;
th.C = nML.C; 
th.Omega = nML.Omega;
th.D = nML.D;

deter = nML.deter;

[n,s]=size(th.K);

par = th2param(th,c,1);
parom = extr_lowtri(th.Omega);

restrict.det_res=0;

[th] = param2th(par,size(th.K,2),n,c,restrict);

th.D = nML.D;
switch deter{1}
    case 'AR'
        th.B = th.K*th.D;
end

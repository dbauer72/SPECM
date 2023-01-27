function nis = poly2ss_th(ni);
% poly2ss_th transforms a system in theta format into a state space system.
%
% SYNTAX: th = poly2ss_th(th);
%
% INPUT: th ... theta structure
%
% OUTPUT: th ... theta structure.
%
% REMARK: converts only the stochastic system (a,b) into state space form.
%
% AUTHOR: dbauer, 14.1.2020

if (strcmpi(ni.which,'SS'))
   disp('Already in state space form!');
   nis = ni;
   return;
end;

% ---- write companion form ----%
a = ni.a;
b = ni.b; 
d= ni.d; 

s = size(a,1); % dimension of endogenous vars. 

if rem(size(b,2),s)
    disp('poly2ss: Wrong dimension of b(z)!');
    return;
end;
if rem(size(a,2),s)
    disp('poly2ss: Wrong dimension of a(z)!');
    return;
end;

p = round(size(a,2)/s-1);
q = round(size(b,2)/s-1);

r = 0; 
m = ni.m;
if ~isempty(d)
    r = round(size(d,2)/m-1);
end; 

if (r>0)    
    n = (p+q)*s+r*m;
    if (q>0)
        A = [-a(:,s+1:end),d(:,m+1:end),b(:,s+1:end);eye((p-1)*s),zeros((p-1)*s,(1+q)*s+r*m);zeros(r*m,p*s),[zeros(m,r*m);eye((r-1)*m),zeros((r-1)*m,m)],zeros(r*m,q*s);zeros(q*s,p*s+r*m),[zeros(s,q*s);eye(s*(q-1)),zeros(s*(q-1),s)]];
        K = [eye(s);zeros(s*(p-1)+r*m,s);eye(s);zeros(s*(q-1),s)];
    else
        A = [-a(:,s+1:end),d(:,m+1:end);eye((p-1)*s),zeros((p-1)*s,(1+q)*s+r*m);zeros(r*m,p*s),[zeros(m,r*m);eye((r-1)*m),zeros((r-1)*m,m)]];
        K = [eye(s);zeros(s*(p-1)+r*m,s)];
    end
    C = A(1:s,:);
    B = [zeros(p*s,m);eye(m);zeros((r-1)*m+q*s,m)];
    D = d(:,1:m);    
else
    m = 0;
    n = (p+q)*s;
    if (q>0)
        A = [-a(:,s+1:end),b(:,s+1:end);eye((p-1)*s),zeros((p-1)*s,(1+q)*s);zeros(q*s,p*s),[zeros(s,q*s);eye(s*(q-1)),zeros(s*(q-1),s)]];
        K = [eye(s);zeros(s*(p-1),s);eye(s);zeros(s*(q-1),s)];
    else
        A = [-a(:,s+1:end);eye((p-1)*s),zeros((p-1)*s,s)];
        K = [eye(s);zeros(s*(p-1),s)];
    end
    C = A(1:s,:);
    D = zeros(s,0);
    B = zeros(n,0);
end

% obtain Hankel matrix. 
H = myhank(A,[K,B],C,3*n,3*n);
% numerical rank of Hankel matrix. 
% correction, if numerical rank is smaller than analytical. 
n = min(n,rank(H));

nmax = 3*n-2;
p = s+m;
[U,S,V] = svd(H);
nis = theta_urs(); 


if n
	Control = S(1:n,1:n)*V(:,1:n)';
	O = U(:,1:n);
	
	nis.A = O(1:nmax*s,:)\O(s+[1:nmax*s],:);
	nis.K = Control(1:n,1:s);
	nis.C = O(1:s,:);
	nis.B = Control(1:n,s+[1:m]);
	nis.D = zeros(s,m);
	nis.which = 'SS';
    nis.num_param = n*(2*s+m)+s*m; 
else
   nis.A = zeros(0,0);
   nis.B = zeros(0,m);
   nis.C = zeros(s,0);
   nis.D = zeros(s,m);
   nis.K = zeros(0,s);
   nis.which = 'SS';
   nis.num_param = 0; 
end;


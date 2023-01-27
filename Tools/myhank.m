function h = myhank(A,B,C,i,j)
%	creates the Hankel matrix of impulse responses
%	of a state space system of the form:
%		xt+1 = Axt + B[zt;et]
%		yt   = Cxt (+D zt + et)
%
% SYNTAX: function h = myhank(A,B,C,i,j)
%
% INPUTS:   A,B,C ... triple of state space system matrices of order n. 
%                     dimensions: C: sxn, B: nxm. 
%	        i     ... number of block rows, default 3*n
%	        j     ... number of block cols, default 3*n
%
% OUTPUT:  H  ... i*s x j*m  Hankel matrix. 
%
% REMARK:	works for any dimensions of the system matrices!
%
% AUTHOR:	dbauer, 27.9.95

% --------  determine the essential dimensions --------

p = size(C,1);
m = size(B,2);
n = size(A,1); 

if (nargin < 5)
	if (nargin < 4)
		i = 3*n;
	end
	j = i;
end

if (i<0)
	i = 3*n;
end

if (j<0)
	j =3*n;
end

h = zeros(p*i,m*j);

% --------  calculate the impulse response sequence --------

r = zeros(p,(i+j-1)*m);
for k = 0:(i+j-2)
	r(:,(k*m)+1:m*(k+1)) = C*A^(k)*B;
end

% --------  write the impulse response into h ---------

for k = 0:(i-1)
	h((k*p)+1:p*(k+1),:) = r(:,(k*m)+1:(k+j)*m);
end


% function nip = ss2poly(ni);
%
% transforms a state space description into a 
% echelon ARMAX description
%
% dbauer, 30.3.2000

function nip = ss2poly(ni);

if (strcmpi(ni.which,'poly'))
   msgbox('Already in polynomial form!');
   nip = ni;
   return;
end;
[p,n] = size(ni.C);
[n,m] = size(ni.B);

nip = ni;
nip.which = 'poly';

% ----- calculate the Hankel matrix ----
H = myhank(ni.A,[ni.B,ni.K],ni.C,n+1,n);

if (rank(H(1:n,:)) == n)|(p==1)
   if (n>p)
      kind=ones(1,p)*fix(n/p);remnp=rem(n,p);                          

      if remnp>0,kind(1:remnp)=kind(1:remnp)+1;end;
      nn = [];
      for i=1:p
            if (kind(i))
               nn = [nn,i:p:kind(i)*p];
            end;
      end;
   else
      kind = [ones(1,n),zeros(1,p-n)];
      nn = 1:n;
   end;
else
   check = 1;
   while check
      kindc = inputdlg([' System not in the Kronecker canonical form!, please enter indices']);
      kind = num2str(kindc{1});
      if (length(kind) == p)
         for i=1:p
            if (kind(i))
               nn = [nn,i:p:kind(i)*p];
            end;
         end;
      else
         nn = [];
      end;
      check = (rank(H(nn,1:n))==n);
   end;
end;
for i=1:p
   if (kind(i))
      hind(i) = p*(kind(i))+i;
   else
      hind(i) = i;
   end;
end;

atil = H(hind,:)/H(nn,:);
nip.a = [eye(p),zeros(p,(max(kind))*p)];
ab = zeros(p,n*p);
ab(:,hind)=eye(p);
ab(:,nn) = -atil;

n = max(kind);

%for j=1:n
%   L((j-1)*p+[1:p],1:m) = H((n-j)*p+[1:p],1:m);
%   K((j-1)*p+[1:p],[1:p]) = H((n-j)*p+[1:p],m+[1:p]);
%end;

for i=1:n+1
   if m
      bb(:,(i-1)*m+[1:m]) = ab(:,i*p+1:(n+1)*p)*H(1:(n-i+1)*p,1:m);
   end;
   
   cb(:,(i-1)*p+[1:p]) = ab(:,i*p+1:(n+1)*p)*H(1:(n-i+1)*p,m+[1:p]);
end;


cb = cb+ab(:,1:size(cb,2));
for i=1:p
   for j=1:(kind(i)+1)
      nip.a(i,(j-1)*p+[1:p]) = ab(i,(kind(i)+1-j)*p+[1:p]);
      if m
         nip.d(i,(j-1)*m+[1:m]) = bb(i,(kind(i)+1-j)*m+[1:m]);
      end;
      
      nip.b(i,(j-1)*p+[1:p]) = cb(i,(kind(i)+1-j)*p+[1:p]);
   end;
end;

ina = inv(nip.a(:,1:p));
nip.b = ina*nip.b;
nip.a = ina*nip.a;


if m
	nip.d = ina*nip.d;
   for i=1:(size(nip.a,2)/p)
      nip.d(:,(i-1)*m+[1:m]) = nip.d(:,(i-1)*m+[1:m]) + nip.a(:,(i-1)*p+[1:p])*ni.d;
   end;
   
else
   nip.d = zeros(p,0);
end;


% tol = 10^(-10);
% nip.na = ~((abs(nip.a)<tol) | (abs(nip.a - 1)<tol));
% nip.nb = ~((abs(nip.b)<tol));
% nip.nc = ~((abs(nip.c)<tol) | (abs(nip.c - 1)<tol));
% nip.which = 'poly';





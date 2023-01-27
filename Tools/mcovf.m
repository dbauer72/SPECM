function R=mcovf(z,M)
%COVF2   Computes covariance function estimates
%
%	see also R = covf(Z,M)
%
%   works also for NaN's (patch by DBauer,22.7.99)
%
%	


[N,nz]=size(z);
z=[z;zeros(M,nz)];
j=[1:N];
for k=1:M
	if any(isnan(z))
		for l=1:nz
			for i=1:nz
				a(l,i) = nanmean(z(j+k-1,l).*z(j,i));
			end;
		end;
	else
		a=z(j+k-1,:)'*z(j,:)/N;
	end;
	R(:,k)=a(:);
end

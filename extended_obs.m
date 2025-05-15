function ob = extended_obs(a,c,s)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Allocate OB and compute each C A^k term
[ny, nx] = size(c);
ob = zeros(s*ny,nx);
ob(1:ny,:) = c;
for k=1:s-1
  ob(k*ny+1:(k+1)*ny,:) = ob((k-1)*ny+1:k*ny,:) * a;
end
end
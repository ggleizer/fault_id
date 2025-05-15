function X = blkhankel(x,s)
%BLKHANKEL Creates a block Hankel matrix of order s given data x
%   Assumes size(x) = [N,n], and s < N
% 
[N,n] = size(x);
% xc = num2cell(x',[1,n]);  % Highly inefficient!!
% idH = hankel(1:s,s:N);
% X = cell2mat(xc(idH));
xv = x';
xv = xv(:);
idH = (1:s*n)' + (0:n:(N-s)*n);           % Hankel subscripts
X = xv(idH);
end
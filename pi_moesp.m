function [A,B,C,D] = pi_moesp(u, y, s, nx)

% FIXME: there is an implementation error wrt to nu > 1

% Build hankel matrices
Y = blkhankel(y, 2*s);
U = blkhankel(u, 2*s);
ny = size(y,2);
nu = size(u,2);
Nt = size(y,1);

% Past and future
Yf = Y(ny*s+1:end,:);
Up = U(1:nu*s,:);
Uf = U(nu*s+1:end,:);

% LQ decomposition
[~, R] = qr([Uf; Up; Yf]',"econ");
% Q = Q';  % unused
R = R';
R11 = R(1:s*nu,1:s*nu);
%R12 = R(s*nu+1:2*s*nu,1:s*nu);  % unused
%R22 = R(s*nu+1:2*s*nu,s*nu+1:2*s*nu);  % unused
R31 = R(2*s*nu+1:end,1:s*nu);
R32 = R(2*s*nu+1:end,s*nu+1:2*s*nu);

% Recover A,C
[Un, Sn, ~] = svd(R32/sqrt(Nt),"econ");  % TODO: a routine to estimate
% the order of the system from Sn
% stem(diag(Sn));

UN = Un(:,1:nx);
C = UN(1:ny,:);
A = UN(1:end-ny,:)\UN(ny+1:end,:);

% Recover, B, D, x
Unperp = Un(:,nx+1:end);
np = size(Unperp,2);
Xi = (Unperp'*R31)/R11;
% Eq. (45)
XiM = reshape(Xi,nu,s*np)';
Unperptran = Unperp';
UH = zeros(s*np,s*ny);
for i = 1:ny
    uny = [Unperptran(:,i:ny:end), zeros(np,s-1)];
    UNYH = blkhankel(uny',s);
    UH(:,(1:s)+s*(i-1)) = UNYH;
end
idx = reshape(1:s*ny,s,ny)';
idx = idx(:);
UH = UH(:,idx);
IU = blkdiag(eye(ny), UN(1:(s-1)*ny,:));
UHIU = UH*IU;
DB = UHIU\XiM;
D = DB(1:ny,:);
B = DB(ny+1:end,:);
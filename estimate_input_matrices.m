function BD = estimate_input_matrices(Uy,O,ny,nfr,i_begin)
%ESTIMATE_INPUT_MATRICES Finds B,D given A,C and range(Y)
%   Given extended observability matrix O, output dimension ny, and a basis
%   for the column space of Y (Hankel matrix of the output y),
%   estimate_input_matrices computes a matrix BD that spans the space of
%   matrices [B; D] for which there exists input u and initial state x0
%   that generate y.
%   
%   BD = ESTIMATE_INPUT_MATRICES(Uy,O,ny,i_begin) considers the estimation
%   of the bottom submatrix of the Markov Parameter matrix, from i_begin to
%   s. In this case, Uy corresponds to a base of the column space of
%   Y_{i_begin:s}.

nx = size(O,2);
s = size(O,1)/ny;
ranky = size(Uy,2);

if nargin <=3
    nfr = 1;
end
if nargin <= 4
    i_begin = 1;
end

% Basic error checks
assert(size(Uy,1) + (i_begin-1)*ny == size(O,1), "Number of rows of Uy and O must be compatible. Are you sure you have used the same horizon 's' for Y and O?");
assert(round(s)==s, "The number of rows of O must be a multiple of ny.")


% Step 1: Block an LTI Toeplitz factorization
% I.e., find a matrix Q such that
%   UQ = [D         0         0 0 ... 0;
%        [CB        D         0 0 ... 0;
%        [CAB       CB        D 0 ... 0;   <--- i_begin
%            :
%        [CA^(s-2)B CA^(s-3)B     ... D];

% Note that
%   UQ = [U1Q1      U1Q2      U1Q3  ... ;
%        [U2Q1      U2Q2      U2Q3  ... ;
%        [U3Q1      U3Q2      U3Q3  ... ;  <--- i_begin
%            :
%        [UsQ1      UsQ2      ...   UsQs];
% where U is divided in block rows and Q in block columns of conformal
% dimension ny.

% Step 1: Find a Toeplitz factorization with zero upper block triangle,
% i.e., for all i = i_begin..s-1 and j = 1..s we have 
%       U_iQ_j = 0,                if j > i   ( (s-i_begin+1)*(s-i_begin)/2 cases)
%       U_iQ_j = U_{i+1}Q_{j+1},   otherwise  ( s*(s-1)/2 - (i_begin)(i_begin-1)/2  cases)
%       O(1:s-1)*P = Uy(2:s)*Q_1  <- Note that P gives the range of B.
% The last condition is for compatibility with A,C: Uy*Q1 = [D; CB; CAB; ...]

% Note: when i_begin > 1, we need additional conditions to ensure
% compatibility with A, C. E.g., s = 5, i_begin = 4; then
%   UQ = [CA2B      CAB       CB       D       0;
%        [CA3B      CA2B      CAB      CB      D]
%   With the previous constraints, we would only ensure that the blocks
%   11, 21 and 22 are compatible with A,C, but not the rest. The most
%   general (while still economic) approach is to add constraints relative
%   to the last row, not the first column. With this, Toeplitz ensures that
%   the whole matrix is an MP matrix compatible with A,C. These conditions
%   are, by equating with
%   UQ = [U1Q1      U1Q2      ...          ;
%                     :
%        [UsQ1      UsQ2      ...      UsQs]
%   for j = 1:s-1
%       O(s-j)*P = U_{s}Q_{j}

% There are s*(s-i_begin) block constraints for Toeplitz, hence we build a
% matrix M with s*(s-i_begin)*ny rows and s*ranky columns. This gives MQ =
% 0. Adding LTI with observability matrix O, we have
%   [M                    0][Q] = 0
%   [-blkdiag(U)   0 O_flip][P]

% See paper: ranky = nx + s*nu - nz >= s*nu for sufficiently large s. So
% both are in the orther of s^2. However, on a system with more outputs
% than inputs, as s grows, the matrix becomes tall and has no structural
% null space. Apparently we want s to be large enough such that this
% happens.

M = zeros((s*(s-i_begin) + s-1)*ny, s*ranky + nx);
% Toeplitz "constraints"
k = 1;  % constraint counter
for i = i_begin:s-1
    iu = i - i_begin + 1;   % Uy is the i_begin:s tail of the original Uy.
    for j = 1:s
        if j > i
            M((k-1)*ny+1:k*ny, (j-1)*ranky+1:j*ranky) = Uy((iu-1)*ny+1:iu*ny,:);
        else
            M((k-1)*ny+1:k*ny, (j-1)*ranky+1:j*ranky) = Uy((iu-1)*ny+1:iu*ny,:);
            M((k-1)*ny+1:k*ny, j*ranky+1:(j+1)*ranky) = -Uy(iu*ny+1:(iu+1)*ny,:);
        end
        k = k+1;
    end
end
% LTI constraint
M((k-1)*ny+1:end, 1:ranky*(s-1)) = -kron(eye(s-1),Uy(end-ny+1:end,:));
blkflipud_idx = floor(((s-1)-1/ny):-1/ny:0)*ny + repmat(1:ny,1,s-1);
Oflipud = O(blkflipud_idx,:);
M((k-1)*ny+1:end, end-nx+1:end) = Oflipud;

% Old LTI constraint
% M((k-1)*ny+1:end, 1:ranky) = -Uy(ny+1:end,:);
% M((k-1)*ny+1:end, end-nx+1:end) = O(1:end-ny,:);

Q = null(M);
if isempty(Q)
    [~,res,Q] = svds(M,nfr,'smallest');  % FIXME: replace 1 by nu (estimated using ranks)
    %warning('There is no Toeplitz-(C,A) decomposition of Uy. Chosing the best approximation with residual error %g', trace(res))
end

Qs = Q(end-ranky-nx+1:end-nx,:);  % We only need Qs to proceed: Us*Qs = D;
Br = Q(end-nx+1:end,:);  % By construction, B in range(P)

% Step 2. Compute B and D
Dr = Uy(end-ny+1:end,:)*Qs;
BDr = [Br; Dr];
% Reduce/orthogonalize basis
BD = orth(BDr);

end
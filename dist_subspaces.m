function d = dist_subspaces(A, B)
% Computes the distance between subspaces spanned by A and B using
% principal angles as per Ye and Lim (2016)

% Ye, K., & Lim, L. H. (2016). Schubert varieties and distances between
%     subspaces of different dimensions. SIAM Journal on Matrix Analysis
%     and Applications, 37(3), 1176-1197.

% Theorem 7 (iii)
A = orth(A);
B = orth(B);

k = size(A,2);
l = size(B,2);

s = svd(A'*B);
s = s(1:min(k,l));
theta = acos(min(s,1));  % Prevent values slightly bigger than 1 for numerical stability

d = norm(theta);
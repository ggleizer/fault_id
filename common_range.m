function X = common_range(varargin)

N = length(varargin);
n = size(varargin{1}, 1);
orths = cellfun(@(x) orth(x), varargin, 'UniformOutput', false);
ncols = cellfun(@(x) size(x,2), orths);

M = BlockMatrix(n*ones(1, N-1), ncols);
for i = 1:N-1
    M.setBlock(i,i, orths{i});
    M.setBlock(i,i+1, -orths{i+1});
end

V = null(M.fullMatrix, 1e-10);
X = M(1,1)*V(1:ncols(1), :);
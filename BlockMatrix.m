classdef BlockMatrix < handle
    properties
        fullMatrix       % The full matrix
        nrows            % Vector of row block sizes
        ncols            % Vector of column block sizes
    end

    methods
        % Constructor
        function obj = BlockMatrix(varargin)
            if nargin == 2
                % Case 2: Only block structure provided (nrows, ncols)
                obj.nrows = varargin{1};
                obj.ncols = varargin{2};
                % Create a zero matrix with the specified block structure
                obj.fullMatrix = zeros(sum(obj.nrows), sum(obj.ncols));
            elseif nargin == 3
                % Case 1: Full matrix and block structure provided
                obj.fullMatrix = varargin{1};
                obj.nrows = varargin{2};
                obj.ncols = varargin{3};
            else
                error('Invalid number of input arguments. Use BlockMatrix(matrix, nrows, ncols) or BlockMatrix(nrows, ncols).');
            end
        end

        % Method to get a block (i,j)
        function block = getBlock(obj, i, j)
            % Calculate row and column indices for the block
            rowStart = sum(obj.nrows(1:i-1)) + 1;
            rowEnd = rowStart + obj.nrows(i) - 1;
            colStart = sum(obj.ncols(1:j-1)) + 1;
            colEnd = colStart + obj.ncols(j) - 1;

            % Extract the block
            block = obj.fullMatrix(rowStart:rowEnd, colStart:colEnd);
        end

        % Method to set a block (i,j)
        function setBlock(obj, i, j, block)
            % Calculate row and column indices for the block
            rowStart = sum(obj.nrows(1:i-1)) + 1;
            rowEnd = rowStart + obj.nrows(i) - 1;
            colStart = sum(obj.ncols(1:j-1)) + 1;
            colEnd = colStart + obj.ncols(j) - 1;

            % Set the block
            obj.fullMatrix(rowStart:rowEnd, colStart:colEnd) = block;
        end

        % Override subsref to handle block access, setting, and slicing
        function varargout = subsref(obj, S)
            switch S(1).type
                case '()'
                    % Handle block access or slicing
                    if numel(S(1).subs) == 2
                        % Extract row and column indices
                        rowIdx = S(1).subs{1};
                        colIdx = S(1).subs{2};

                        % Convert block indices to matrix indices
                        if isscalar(rowIdx) && isscalar(colIdx)
                            % Single block access
                            varargout{1} = obj.getBlock(rowIdx, colIdx);
                        else
                            % Slicing: convert block indices to matrix indices
                            rowStart = sum(obj.nrows(1:rowIdx(1)-1)) + 1;
                            rowEnd = sum(obj.nrows(1:rowIdx(end)));
                            colStart = sum(obj.ncols(1:colIdx(1)-1)) + 1;
                            colEnd = sum(obj.ncols(1:colIdx(end)));

                            % Extract the submatrix
                            varargout{1} = obj.fullMatrix(rowStart:rowEnd, colStart:colEnd);
                        end
                    else
                        error('Invalid indexing. Use (i,j) or (i1:i2, j1:j2).');
                    end
                otherwise
                    % Fallback to default behavior
                    [varargout{1:nargout}] = builtin('subsref', obj, S);
            end
        end

        % Override subsasgn to handle setting blocks
        function obj = subsasgn(obj, S, value)
            switch S(1).type
                case '()'
                    % Handle block assignment
                    if numel(S(1).subs) == 2
                        % Extract row and column indices
                        rowIdx = S(1).subs{1};
                        colIdx = S(1).subs{2};

                        if isscalar(rowIdx) && isscalar(colIdx)
                            % Single block assignment
                            obj.setBlock(rowIdx, colIdx, value);
                        else
                            % Slicing: convert block indices to matrix indices
                            rowStart = sum(obj.nrows(1:rowIdx(1)-1)) + 1;
                            rowEnd = sum(obj.nrows(1:rowIdx(end)));
                            colStart = sum(obj.ncols(1:colIdx(1)-1)) + 1;
                            colEnd = sum(obj.ncols(1:colIdx(end)));

                            % Assign the submatrix
                            obj.fullMatrix(rowStart:rowEnd, colStart:colEnd) = value;
                        end
                    else
                        error('Invalid indexing. Use (i,j) or (i1:i2, j1:j2).');
                    end
                otherwise
                    % Fallback to default behavior
                    obj = builtin('subsasgn', obj, S, value);
            end
        end

        % Method to sum two block matrices
        function result = plus(obj1, obj2)
            % Check if block structures are the same
            if ~isequal(obj1.nrows, obj2.nrows) || ~isequal(obj1.ncols, obj2.ncols)
                error('Block structures must be the same for addition.');
            end
            result = BlockMatrix(obj1.fullMatrix + obj2.fullMatrix, obj1.nrows, obj1.ncols);
        end

        % Method to multiply two block matrices
        function result = mtimes(obj1, obj2)
            % Check if block structures are conformal
            if ~isequal(obj1.ncols, obj2.nrows)
                error('Block structures are not conformal for multiplication.');
            end
            result = BlockMatrix(obj1.fullMatrix * obj2.fullMatrix, obj1.nrows, obj2.ncols);
        end

        function disp(obj)
            % Display the full matrix with block structure
            [rows, cols] = size(obj.fullMatrix);
            rowBreaks = cumsum(obj.nrows(1:end-1)); % Exclude the last break
            colBreaks = cumsum(obj.ncols(1:end-1)); % Exclude the last break

            % Define fixed width for each cell
            cellWidth = 8; % Adjust this value as needed

            % Display the matrix with separators
            for i = 1:rows
                % Add horizontal line after each block row (except the last)
                if ismember(i, rowBreaks + 1) % Shift by 1 to place after the block
                    disp(repmat('-', 1, cols * (cellWidth + 1))); % Add horizontal line
                end
                for j = 1:cols
                    % Add vertical line after each block column (except the last)
                    if ismember(j, colBreaks + 1) % Shift by 1 to place after the block
                        fprintf('|'); % Add vertical line
                    end
                    % Print the matrix element with fixed width
                    fprintf(sprintf('%%%d.4f ', cellWidth), obj.fullMatrix(i, j));
                end
                fprintf('\n');
            end
        end
    end
end
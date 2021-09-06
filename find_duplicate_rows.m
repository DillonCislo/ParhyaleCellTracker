function [ dupl, C ] = find_duplicate_rows( A, setOrder )
%FIND_DUPLICATE_ROWS Finds duplicate rows in numeric matrices
%
%   INPUT PARAMETERS:
%
%       - A:            A numeric array
%
%       - setOrder:     Returns the unique rows of A with a specified order
%                       - {'sorted'} Ascending order by columns
%                       - 'stable' Stable ordering
%
%   OUTPUT PARAMETERS:
%
%       - dupl:         A struct with two fields and N elements (one for
%                       each one of the N duplicate rows)
%                       - 'val': The value of the duplicated row
%                       - 'idx': The row IDs of all those rows sharing the
%                                value in 'val'
%
%       - C:            The unique rows of A
%                           
% by Dillon Cislo 01/22/2020

if (nargin < 2), setOrder = 'sorted'; end
if ~(strcmp(setOrder, 'sorted') || strcmp(setOrder, 'stable'))
    error('Invalid set order supplied for unique row output');
end

[C, ia, ~] = unique( A, 'rows', setOrder );

if size(A,1) == size(C,1)
    % disp('There are no duplicate rows!');
    dupl = {};
    return;
end

rep_idx = setdiff(1:size(A,1), ia);
rep_val = unique( A(rep_idx,:), 'rows');

dupl_val = cell( size(rep_val,1), 1 );
dupl_idx = cell( size(rep_val,1), 1 );

for i = 1:size(rep_val,1)
    
    dupl_val{i} = rep_val(i,:);
    
    diffA = A - repmat( rep_val(i,:), size(A,1), 1 );
    dupl_idx{i} = find( abs(sum( diffA, 2 )) < 10*eps );
    
end

dupl = struct( 'val', dupl_val, 'idx', dupl_idx );

end


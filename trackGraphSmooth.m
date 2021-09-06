function smoothFields = trackGraphSmooth(G, fieldNames, varargin)
%TRACKGRAPHSMOOTH This function acts a type of discrete "Gaussian"
%smoothing in time on a set of fields from an input tracking graph.
%
%   NOTE: FOR NOW IT IS IMPLICITLY ASSUMED THAT THE TIME STEP BETWEEN ALL
%   NODES IS IDENTICAL
%
%   INPUT PARAMETERS:
%
%       - G:            A cell tracking digraph.  Must contain a field 'T'
%                       corresponding to the time point of each node.  It
%                       is assumed that each node has one and only one
%                       predecessor, but each node may have multiple
%                       successors (i.e., due to cell division)
%
%       - fieldNames:   N-element cell array of character vectors
%                       containing the names of fields to smooth
%
%       (Name, Value)-Pair Arguments:
%           
%               - 'NumIter':        The number of smoothing iterations {10}
%
%               - 'NumNeighbors':   The number of neighbors used to smooth
%                                   the fields.  Must be an even number
%                                   (i.e. NN/2 predecessors and NN/2
%                                   successors). Weighting of each neighbor
%                                   will be given according to the
%                                   corresponding row of Pascal's triangle
%                                   {4}
%
%               - 'CensorNodes':    A logical vector or list of node IDs
%                                   that determine which nodes will be
%                                   censored from consideration during the
%                                   smoothing operation {[]}
%
%               - 'CensorDiv':      If true, all nodes in a temporal window
%                                   of specified length around division
%                                   events will be censored. {'false'}
%
%               - 'DivWindow':      The size of the temporal window
%                                   around division events within which
%                                   nodes will be censored if desired {6}
%
%               - 'CensorTimes':    A list of time points to censor {[]}
%
%   OUTPUT ARGUMENTS:
%
%       - smoothFields: Nx1 cell array. Each entry is the smoothed fields
%                       associated to the corresponding entry in
%                       'fieldNames'
%
% by Dillon Cislo 02/07/2020

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

if (nargin < 1), error('Please supply tracking graph'); end
if (nargin < 2), error('Please supply field names to smooth'); end

% Validate the Input Tracking Graph ---------------------------------------
assert( isa(G, 'digraph'), ...
    'Cell tracking must be supplied as a digraph object' );

% Make sure that no nodes have multiple predecessors
oneParent = isequal( unique( G.Edges.EndNodes(:,2), 'stable' ), ...
    G.Edges.EndNodes(:,2) );
assert( oneParent, ...
    'Nodes in the supplied tracking graph have multiple parents' )

propNames = G.Nodes.Properties.VariableNames;
assert( ismember( 'T', propNames ), ...
    'Tracking graph must have a "T" field denoting time' );

% Validate the Fields to be Smoothed --------------------------------------
assert( all( cellfun( @(x) isa(x, 'char'), fieldNames ) ), ...
    [ 'Fields to be smoothed must be supplied as a cell' ...
    'array of character vectors' ] );

assert( ~ismember('T', fieldNames), 'The time field cannot be smoothed' );
assert( all(ismember(fieldNames, propNames)), 'Invalid field names' );

% Extract the input fields
inputFields = cell(size(fieldNames));
for i = 1:numel(inputFields)
    inputFields{i} = G.Nodes.(fieldNames{i});
    assert( isnumeric(inputFields{i}), 'Only numeric can be smoothed' );
    assert( ismatrix(inputFields{i}), ...
        'This function only supports matrix fields (for now)' ); 
end

% Process Optional Inputs -------------------------------------------------
numIter = 10;
numNeighbors = 4;
censorIDx = [];
censorDiv = false;
divWindow = 5;
censorTimes = [];

for i = 1:length(varargin)
    if isa(varargin{i},'double') 
        continue;
    end
    if isa(varargin{i},'logical')
        continue;
    end
    if ~isempty(regexp(varargin{i},'^[Nn]um[Ii]ter','match'))
        numIter = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Nn]um[Nn]eighbors','match'))
        numNeighbors = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Cc]ensor[Nn]odes','match'))
        censorIDx = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Cc]ensor[Dd]iv','match'))
        censorDiv = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Dd]iv[Ww]indow','match'))
        divWindow = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Cc]ensor[Tt]imes','match'))
        censorTimes = varargin{i+1};
    end
end

validateattributes( numIter, {'numeric'}, ...
    {'scalar', 'finite', 'real', 'integer', 'nonnegative'} );

validateattributes( numNeighbors, {'numeric'}, ...
    {'scalar', 'finite', 'real', 'integer', 'nonnegative', 'even'} );

validateattributes( censorDiv, {'logical'}, {'scalar'} );

validateattributes( divWindow, {'numeric'}, ...
    {'scalar', 'finite', 'real', 'integer', 'nonnegative', 'even'} );

if ~isempty(censorTimes)
    validateattributes( censorTimes, {'numeric'}, ...
        {'vector', 'finite', 'real', 'integer', 'nonnegative'});
end

if any(~ismember(censorTimes, unique(G.Nodes.T)))
    error('Invalid list of censored time points');
end

censorIDx = censorIDx(:); % Convert to column vector
if islogical(censorIDx)
    
    if (numel(censorIDx) ~= numnodes(G))
        error('Logical vector of censored nodes is improperly sized');
    end
    
    % Convert to list of node IDs
    censorIDx = find(censorIDx);
    
else
    
    validateattributes( censorIDx, {'numeric'}, ...
        {'vector', 'finite', 'real', 'integer', ...
        'positive', '<=', numnodes(G)} );
    
    % Remove duplicate entries
    censorIDx = unique(censorIDx);
    
end

% Exit the function in the trivial cases
if ( (numIter == 0) || (numNeighbors == 0) )
    smoothFields = inputFields;
    return;
end

%--------------------------------------------------------------------------
% HANDLE DIVISION CENSORING
%--------------------------------------------------------------------------
% NOTE: This current implementation is inefficient and may be slow if
% either the number of division events or the size of the division
% censoring window is large. If necessary, it should be possible to simply
% incorporate it into the neighbor determination process with minimal
% overhead

sources = G.Edges.EndNodes(:,1); % The source nodes of each edge
sinks = G.Edges.EndNodes(:,2); % The sink nodes of each edge

% Determine the edges whose sources have only a Single Successor
% i.e., edges corresponding to a division event
SS = findUniqueEntryLocations( sources );

if censorDiv
    
    % Censor the time point leading up to a division event ----------------
    parentIDx = unique(sources(~SS));
    parentIDx = parentIDx(:);

    censorParentIDx = zeros(numel(parentIDx), divWindow/2);
    censorParentIDx(:,1) = parentIDx;
    
    for i = 1:numel(parentIDx)
        for j = 2:(divWindow/2)
            
            predID = predecessors(G, censorParentIDx(i,j-1));
            if ~isempty(predID)
                censorParentIDx(i,j) = predID;
            end
            
        end
    end
    
    censorParentIDx(censorParentIDx == 0) = [];
    censorParentIDx = censorParentIDx(:);
    
    % Censor the time points following a division event -------------------
    % Perform a depth-first search to find the successors of a division

    % The search stack
    childIDx = unique(sinks(~SS));
    childIDx = childIDx(:); 
    
    % A logical vector indicating if nodes have been visited
    visitedChild = false(numnodes(G), 1);
    
    % A vector specifying the depth of nodes relative to the child nodes
    childDepth = zeros(numnodes(G), 1);
    
    % The censored ID vector
    censorChildIDx = [];
    
    % The DFS search
    while ~isempty(childIDx)
        
        % Pop the current child ID off of the search stack
        childID = childIDx(1);
        childIDx(1) = [];
        
        % Check if the node has already been visited
        if ~visitedChild(childID)
            
            % Label the node as discovered
            visitedChild(childID) = true;
            
            % Check if the node has sufficiently shallow depth
            if (childDepth(childID) <= (divWindow/2))
                
                % Find the successors of the node
                succID = successors(G, childID);
                succID = succID(:);
                
                % Update the depth of the successors
                childDepth(succID) = childDepth(childID) + 1;
                
                % Push the successors on top of the search stack
                childIDx = [ succID; childIDx ];
                
                % Add the current child ID to the censored ID vector
                censorChildIDx = [ censorChildIDx; childID ];
                
            end
            
        end
        
    end
    
    % Concatenate censored times to the extant censored nodes list --------
    censorIDx = [ censorIDx; censorParentIDx; censorChildIDx ];

end

% Censor Specified Time Points --------------------------------------------
if ~isempty(censorTimes)
    for i = 1:numel(censorTimes)
        
        nodeIDx = find(G.Nodes.T == censorTimes(i));
        censorIDx = [ censorIDx; nodeIDx(:) ];
        
    end
end

% Remove repeat entries
censorIDx = unique(censorIDx);
censorIDx = censorIDx(:);

%--------------------------------------------------------------------------
% COMPUTE THE NEIGHBORS OVER WHICH TO SMOOTH
%--------------------------------------------------------------------------
% For the purposes of smoothing the branches will terminate at division
% events in both the forward and backward directions. Censored points will
% not cause a branch to terminate but will be removed from consideration
% later

% A array to hold the IDs of the neighbors
% [ -NN/2 (-NN/2+1) ... -1 NODE 1 ... (NN/2-1) NN/2 ]
neighborIDs = nan( numnodes(G), numNeighbors+1 );
neighborIDs( :, (numNeighbors/2)+1 ) = (1:numnodes(G)).';

% Search for Successors ---------------------------------------------------

% Assign the first successor
neighborIDs( sources(SS), (numNeighbors/2)+2 ) = sinks(SS);
firstSucc = neighborIDs(:, (numNeighbors/2)+2 );

% Assign all higher order successors
for i = (numNeighbors/2+3):(numNeighbors+1)
    
    curSource = neighborIDs(:, i-1);
    
    % A logical vector indicating active sources
    extantSourceLoc = ~isnan(curSource);
    
    % The actual node IDs of active sources
    extantSourceNID = curSource(extantSourceLoc);
    
    curSink = nan(size(curSource));
    curSink(extantSourceLoc) = firstSucc(extantSourceNID);
    
    neighborIDs(:, i) = curSink;
    
end

% Search for Predecessors -------------------------------------------------

% Assign the first predecessor
neighborIDs( sinks(SS), numNeighbors/2 ) = sources(SS);
firstPred = neighborIDs( :, numNeighbors/2 );

% Assign all higher order predecessors
for i = (numNeighbors/2-1):-1:1
    
    curSink = neighborIDs(:, i+1);
    
    % A logical vector indicating active sinks
    extantSinkLoc = ~isnan(curSink);
    
    % The actual node IDs of active sinks
    extantSinkNID = curSink(extantSinkLoc);
    
    curSource = nan(size(curSink));
    curSource(extantSinkLoc) = firstPred(extantSinkNID);
    
    neighborIDs(:, i) = curSource;
    
end

% Set any nensored nodes as NaN entries
neighborIDs( ismember(neighborIDs, censorIDx) ) = NaN;

% Restore the ID for the basic node list
neighborIDs( :, (numNeighbors/2)+1 ) = (1:numnodes(G)).';

%--------------------------------------------------------------------------
% CONSTRUCT WEIGHTS FOR SMOOTHING VIA WEIGHTED AVERAGE
%--------------------------------------------------------------------------

W = computePascalRow(numNeighbors);
W = repmat( W, numnodes(G), 1 );

W( isnan(neighborIDs) ) = 0;
WSum = sum(W, 2);

%--------------------------------------------------------------------------
% PERFORM ITERATIVE SMOOTHING
%--------------------------------------------------------------------------

% Remove NaNs for vectorized indexing (these values will not be included in
% the weighted average anyway)
neighborIDs( isnan(neighborIDs) ) = 1;

N = numNeighbors+1;

smoothFields = inputFields;

for i = 1:numIter
    for j = 1:numel(smoothFields) 
        
        curField = smoothFields{j};
        
        if isvector(curField)
            
            smoothField = sum( W .* curField( neighborIDs ), 2 ) ./ WSum;
            
        else
            
            numCol = size(curField, 2);
            
            % smoothField = zeros( numnodes(G), numCol );
            % 
            % for k = 1:numCol
            %     
            %     curCol = smoothField(:,k);
            %     
            %     smoothField(:,k) = ...
            %         sum( W .* curCol( neighborIDs ), 2 ) ./ 2;
            %     
            % end
            
            % OLD METHOD: I really don't get what I was thinking here
            % -------------------------------------------------------
            curField = curField( neighborIDs.', : );
            curField = repmat( reshape(W.', numel(W), 1), ...
                1, size(curField, 2) ) .* curField;
            
            smoothField = zeros( numnodes(G), numCol );
            
            for k = 1:numnodes(G)
               
                smoothField(k,:) = sum( ...
                    curField( (1+N*(k-1)):(N*k), : ), 1 ) ./ WSum(k);
                
            end
            
        end
        
        smoothFields{j} = smoothField;
  
    end
end


end

function pascalRow = computePascalRow(n)
%COMPUTEPASCALROW Computes the nth row of Pascals triangle
%
%   INPUT ARGUMENTS:
%
%       - n:            A non-negative integer
%
%   OUTPUT ARGUMENTS:
%
%       - pascalRow:    The nth row of Pascal's triangle.  Counting starts
%                       at n = 0

if (n == 0), pascalRow = 1; return; end
if (n == 1), pascalRow = [1 1]; return; end

pascalRow = zeros(1, n+1);
pascalRow([1 2]) = [1 1];
for i = 3:(n+1)
    pascalRow(1:i) = ...
        [ 1, ( pascalRow(1:(i-2)) + pascalRow(2:(i-1)) ),  1 ];
end

end

function loc = findUniqueEntryLocations( V )
%FINDUNIQUEENTRYLOCATIONS Find the locations of the unique entries of a
%vector.
%
%   INPUT ARGUMENTS:
%
%       - V:        N-element numeric vector
%
%   OUTPUT ARGUMENTS:
%
%       - loc:      Logical vector the same size as V. Entries are true if
%                   the corresponding entry in V is unique

[ ~, iv, ~ ] = unique( V, 'stable' );

duplVal = V;
duplVal(iv) = [];

loc = ~ismember(V, duplVal);

end



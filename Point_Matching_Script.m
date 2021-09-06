%% Parhyale Cell Tracking =================================================
%
% This is a script to perform automated cell tracking via simple point
% matching. It also includes a section to call a GUI for manual correction
% of the automatic results. The input to this method is a set of images
% (one for each time point), corresponding to your data, and a cell array
% holding a list of unordered 2D cell centers at each time point. The list
% of cell centers may contain an unknown number of erroneous duplicate
% cells at each time point (the potential product of a faulty segmentation
% pipeline). The output is a complete cell tracking stored as a MATLAB
% digraph
%
% by Dillon Cislo 09/05/2021
%
%==========================================================================

clear; close all; clc;

% Add the manual tracking code to the path
addpath(genpath('PATH_TO_PCT/ParhyaleCellTracking/Manual_Tracking_GUI'));

% Define a working directory
[projectDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(projectDir);

% The directory holding the pullback files
% *** INSERT YOUR OWN FILES HERE ***
MIPDir = 'PATH_TO_IMAGE_FILES';

% The naming convection of the pullback files
MIPFile = fullfile( MIPDir, 'cmp_1_1_T%04d.tif' );

% Load the cell center tracks
% This is a Tx1 cell array
% *** INSERT YOUR OWN FILES HERE ***
load('PATH_TO_CELL_CENTERS/CELL_CENTERS.mat');

%%
clear; close all; clc;

% Add the manual tracking code to the path
addpath(genpath('/data/code/ParhyaleCellTracker'));

% Define a working directory
[projectDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(projectDir);

% The directory holding the pullback files
MIPDir = [ '/mnt/crunch/djcislo/MATLAB/Parhyale_Imsane/' ...
    'Parhyale_Imsane_Normal_2019111/' ...
    'Parhyale_SOI_T1_T434_No_Shift_20191107/' ...
    'fields/data_MIP/xy_index/xy' ];

% The naming convection of the pullback files
MIPFile = fullfile( MIPDir, 'cmp_1_1_T%04d.tif' );

% Load the cell center tracks
load(['/mnt/crunch/djcislo/MATLAB/CellTracking/' ...
    'Parhyale_Cell_Tracking_20200121/cell_centers.mat']);

%% ************************************************************************
% *************************************************************************
% *************************************************************************
%           PERFORM AUTOMATED TRACKING VIA POINT MATCHING
% *************************************************************************
% *************************************************************************
% *************************************************************************

%% Add Cells to Formatted Digraph Object ==================================

% A list of the coordinates of all cells from all time points
% NOTE: We have to do this stupid method since there are unknown numbers of
% duplicates at each time point
allCellCenters = [];

% A list of the time points corresponding to each cell
allCellTimes = [];

count = 1; % An indexing variable
for t = 1:numel(cell_centers)
    
    % The number of cells at the current time point
    numCellsT = size(cell_centers{t}, 1);
    
    % Remove duplicate cell entries
    [~, uniqueCC] = find_duplicate_rows( cell_centers{t}, 'stable' );
    
    if size(uniqueCC, 1) ~= numCellsT
        numDupl = numCellsT - size(uniqueCC, 1);
        numCellsT = size(uniqueCC, 1);
        warning( '%d duplicate cells found at entry %d\n', numDupl, t );
    end
    
    % Update the global cell center list
    allCellCenters = [allCellCenters; uniqueCC];
    
    % Update the times in the global list
    allCellTimes =  [allCellTimes; repmat(t, numCellsT, 1)];
    
    % Update the indexing variable
    count = count+numCellsT;
    
end

% Update the time points to reflect that the count may actually starts at a
% timepoint other than 1. In this example, the count starts at 201;
timeShift = 200;
allCellTimes = allCellTimes + timeShift;

clear count t numCellsT uniqueCC

% Generate the Segment Label/Generation fields for the graph
allSegmentLabels = repmat({'none'}, numel(allCellTimes), 1);
allGenerations = ones(numel(allCellTimes), 1);

% A table containing node properties with which to construct the digraph
nodeTable = table( allCellCenters, allCellTimes, ...
    allSegmentLabels, allGenerations, ...
    'VariableNames', {'UPix', 'T', 'Segment', 'Generation'} );


% Construct the graph and update node properties
G = digraph;
G = G.addnode(nodeTable);

clear allCellCenters allCellTimes nodeTable cell_centers
clear allSegmentLabels allGenerations

%% Perform Point Matching =================================================

% All of the time points contained in the segmentation data set
timePoints = 201:434;

for tidx = 2:numel(timePoints)
    
    t = timePoints(tidx);
    fprintf('Processing time point T = %d\n', t);
    
   % Find all of the cell centers at the current time
   curID = find(G.Nodes.T == t);
   curXY = G.Nodes(curID,:).UPix;
   
   % Find all of the cell centers at the previous time
   prevID = find(G.Nodes.T == timePoints(tidx-1));
   prevXY = G.Nodes(prevID,:).UPix;
   
   % Find the ID of the nearest cell at the previous time to each cell
   % at the current time
   parentID = zeros(size(curID));
   
   for i = 1:numel(curID)
       
       nearID = prevXY - repmat( curXY(i,:), numel(prevID), 1 );
       nearID = sqrt( sum( nearID.^2, 2 ) );
       [~, nearID] = min(nearID);
       
       parentID(i) = nearID;
       
   end
   
   % Map cell IDs to row IDs in the digraph node table
   parentID = prevID(parentID);
   
   % Update the edge list in the digraph
   G = addedge(G, parentID, curID);
    
end

clear curID curXY prevID prevXY parentID nearID

%% Perform Manual Tracking ================================================

% Open manual tracking gui
[Gout, divStruct] = parhyale_master_gui(G, MIPFile);


%% ************************************************************************
% *************************************************************************
% *************************************************************************
%                CREATE MOVIE FOR LINEAGE VISUALIZATION
% *************************************************************************
% *************************************************************************
% *************************************************************************

%% Assign Basic Lineage IDs ===============================================

%--------------------------------------------------------------------------
% Determine the 'progenitor' cells which have no parents in the tracking
% graph structure
%--------------------------------------------------------------------------

% Find the in-degree of each node
inDeg = indegree(G);

assert(isequal(unique(inDeg), [0; 1]), 'Invalid tracking structure');

progCells = find(inDeg == 0);

clear inDeg

%--------------------------------------------------------------------------
% Add a lineage ID field to each node of the tracking graph
%--------------------------------------------------------------------------

lineageID = -ones(size(G.Nodes,1),1);
lineageID(progCells) = 1:numel(progCells);

G.Nodes.LineageID = lineageID;

clear lineageID

%--------------------------------------------------------------------------
% Update the lineage ID field of each node using a depth-first search
%--------------------------------------------------------------------------

for i = 1:numel(progCells)
    
    % Get the IDs of all nodes descended from the current progenitor cell
    descNodes = dfsearch(G, progCells(i), ...
        { 'discovernode', 'edgetonew' } );
    descNodes = descNodes.Node;
    descNodes = descNodes( ~isnan(descNodes) );
    
    descNodes( ismember(descNodes, progCells) ) = [];
    
    % Udpate the lineage ID
    G.Nodes(descNodes,:).LineageID = repmat(i, numel(descNodes), 1);
    
end

clear descNodes

assert( ~any(G.Nodes.LineageID < 0), 'Lineage improperly assigned!' );

%% Generate Movie Frames ==================================================

% The naming convection of the pullback files
MIPFile = fullfile( MIPDir, 'cmp_1_1_T%04d.tif' );

% the naming convention of the saved image files
saveFile = [ '/mnt/crunch/djcislo/MATLAB/CellTracking/' ...
    'Parhyale_Cell_Tracking_20200121/Corrected_Tracking_20200213/' ...
    'Corrected_Tracking_T%04d' ];

% A list of colors for each lineage
linColors = distinguishable_colors( numel(progCells), [0 0 0; 1 1 1] );

% A struct to hold the image frames
% F(numel(timePoints)) = struct('cdata', [], 'colormap', [] );

parfor tidx = 1:numel(timePoints)
    
    t = timePoints(tidx);
    fprintf('Now processing timepoints T = %d\n', t);
    
    % Node IDs of the cells at the current time point
    cIDx = find(G.Nodes.T == t);
    
    % The cell centers in pixel space coordinates
    % C2D = [ G.Nodes(cIDx,:).X, G.Nodes(cIDx,:).Y ];
    C2D = G.Nodes(cIDx,:).UPix;
    
    % The lineage IDs of the cells
    LID = G.Nodes(cIDx,:).LineageID;
    
    % Load the current MIP from file
    MIP = imadjust( mat2gray( imread( sprintf( MIPFile, t ) ) ) );
    
    % Capture frames ------------------------------------------------------
    fig = figure('units', 'normalized', 'outerposition', [0 0 1 1], ...
        'visible', 'off');
    
    imshow(MIP)
    set(gca, 'YDir', 'normal');
    hold on
    scatter(C2D(:,1), C2D(:,2), [],  linColors(LID,:), 'filled');
    hold off
    
    % F(tidx) = getframe;
    
    % Write to file -------------------------------------------------------
    imFileName = sprintf( saveFile, t );
    
    export_fig( imFileName, '-tif', '-transparent' );
    
    close(fig);
    
end

disp('DONE');

%% Write to Video File ====================================================

v = VideoWriter('Corrected_Tracking_20200213.avi');
open(v);

for T = 1:numel(timePoints)
    
    writeVideo(v, F(T));
    
end

close(v);

    

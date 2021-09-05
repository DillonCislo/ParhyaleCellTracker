%% PARHYALE_MASTER_GUI_SCRIPT =============================================
%
%   This is pipeline to illustrate the mitotic wave structure in the
%   developing Parhyale germband.  The first primary output of this script
%   is an assignment of cell cycle generation to the segmented cells. Note
%   that this assignment is flawed due to the staggered development of the
%   different parasegment pre-cursor rows and must be manually corrected 
%
%   by Dillon Cislo 06/01/2020
%
%==========================================================================

clear; close all; clc;

% Add some necessary directories to the path ------------------------------
% NOTE: ImSAnE should be added separately using 'setup.m'
% THe ImSAnE installation can be found at /data/code/imsane_for_git 
addpath('/data/code');
addpath('/mnt/crunch/djcislo/MATLAB/CellTracking/');
addpath(genpath('/data/code/euclidean_orbifolds'));
addpath(genpath('/data/code/RicciFlow'));
addpath(genpath('/data/code/TexturePatch'));
addpath(genpath('/data/code/gptoolbox'));
addpath(genpath('/data/code/export_fig'));
addpath(genpath('/data/code/coordinate_transformations'));
addpath(genpath('/data/code/TriStream'));
addpath(genpath('/data/code/ellipse_fitting'));
addpath(genpath('/data/code/DEC'));
addpath(genpath('/data/code/PerceptuallyUniformColormaps'));
addpath(genpath('/data/code/shadedErrorBar'));
addpath(genpath('/data/code/CircStat'));

% Define directories ------------------------------------------------------

% The working directory
[projectDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(projectDir);

% Add the master GUI to the path
addpath(genpath(fullfile(projectDir, 'Parhyale_Master_GUI')));

% Set the MIP file naming convention
MIPFile = ['/mnt/crunch/djcislo/MATLAB/Parhyale_Imsane/' ...
    'Parhyale_Imsane_Normal_2019111/' ...
    'Parhyale_SOI_T1_T434_No_Shift_20191107/fields/data_MIP/xy_index/xy'];

MIPFile = fullfile( MIPFile, 'cmp_1_1_T%04d.tif' );

% The directory to which the output will be saved
saveDir = fullfile( projectDir, ...
    'All_Segments_20200721' );

% Create the output directory if necessary
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

% Load data from file -----------------------------------------------------

% Load tracking graph
% load(fullfile(projectDir, ...
%     'CellTrackingGraph_DynamicNormal_20200823.mat'));

load(fullfile(projectDir, ...
    'CellTrackingGraph_DynamicNormal_RowLabel_T201-385_20200825.mat'));

%%

 [Gout, divStruct] = parhyale_master_gui(G, MIPFile);
 
%% ************************************************************************
% *************************************************************************
%                       VISUALIZE RESULTS
% *************************************************************************
% *************************************************************************

%%

%--------------------------------------------------------------------------
% Set Output File Naming Convention
%--------------------------------------------------------------------------
[saveDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
saveDir = fullfile( saveDir, 'Row_Labels_20200825'  );

if ~exist( saveDir, 'dir' ), mkdir( saveDir ); end

saveFile = fullfile( saveDir, 'Parasegments_T%04d' );

%--------------------------------------------------------------------------
% Extract Parasegment Label/Cell Cycle Generation Lists
%--------------------------------------------------------------------------

segLabels = unique( G.Nodes.Segment );
segLabels( strcmp( segLabels, 'none' ) ) = [];
% segLabels = {'PS2'};

% Determine which nodes should contribute to the coloring scheme
hasLabel = ismember(G.Nodes.Segment, segLabels);

% The number of unique generations
numGen = max(G.Nodes(hasLabel,:).Generation);

% The number of unique labels
numLabels = numel(segLabels);

% Find the nodes corresponding to each label
labelIDx = cell( numLabels, 1 );
for i = 1:numLabels
    labelIDx{i} = find(strcmp( G.Nodes.Segment, segLabels{i} ));
end

%--------------------------------------------------------------------------
% Handle Label/Generation Color Selection
%--------------------------------------------------------------------------
% Each generation up to and including the first round of differential
% cleaveage will be assigned a unique color.  All further rounds of
% differential cleavage will be assigned the same color as the first round

% Assign a unique color to each label
segColors = distinguishable_colors( numLabels + numGen, [0 0 0; 1 1 1] );

genColors = segColors(1:numGen, :);
segColors = segColors((numGen+1):end, :);

% Set erroneous cells to black
% badLabels = {'MB'};
badLabels = {};
badLabel = ismember( segLabels, badLabels);
segColors(badLabel,:) = repmat( [0 0 0], sum(badLabel), 1 );

%--------------------------------------------------------------------------
% Iterate Over Time Points
%--------------------------------------------------------------------------

timePoints = unique(G.Nodes.T);

parfor tidx =  1:numel(timePoints)
    
    t = timePoints(tidx); % The current time point being rendered
    fprintf( 'NOW PROCESSING T = %d\n', t );
    
    %----------------------------------------------------------------------
    % Load/Extract Data
    %----------------------------------------------------------------------
    
    % Load background image from file
    MIP = imadjust( mat2gray( imread( sprintf( MIPFile, t ) ) ) );
    
    % Extract node IDs
    nodeIDx = find(G.Nodes.T == t);
    
    % The cell locations in pixel space
    U = G.Nodes(nodeIDx,:).UPix;
    
    % Extract the generations of the current nodes
    curGen = G.Nodes(nodeIDx,:).Generation;
    
    %----------------------------------------------------------------------
    % Process Voronoi Connectivity
    %----------------------------------------------------------------------
    
    % Construct Voronoi tessellation of the image plane
    delTri = delaunayTriangulation( U );
    [ v, c ] = voronoiDiagram( delTri );
    
     % The cell IDs for each bond
    bondIDx = delTri.edges;
    
    % Determine which cells lie on the boundary of the diagram
    bdyCells = cellfun( @(x) ismember(1, x), c );
    
    % Consider a cell to be a boundary cell if its nth neighbor is also a
    % boundary cell.  Useful for make prettier plots without strange
    % Voronoi cells
    nnBdy = 2;
    for i = 1:nnBdy
        
        bdyCellsTemp = bdyCells;
        for ic = 1:size(c,1)
            
            containsIC = any(bondIDx == ic, 2);
            containsBdy = any(ismember(bondIDx, find(bdyCells)), 2);
            
            if any( containsIC & containsBdy )
                bdyCellsTemp(ic) = true;
            end
        end
        
        bdyCells = bdyCellsTemp;
        
    end
    
    % clear nnBdy bdyCellsTemp ic containsIC containsBdy
    
    maxFaceSize = max(cellfun(@(x) numel(x), c));
    voronoiFace = nan(size(c,1), maxFaceSize);
    for k = 1:size(c,1)
        voronoiFace(k, 1:numel(c{k})) = c{k};
    end
    
    % clear maxFaceSize
    
    % INCLUDE ALL CELLS
    bdyCells = false(size(bdyCells));
    
    %----------------------------------------------------------------------
    % Process Patch Face Colors/Determine Which Nodes Will Be Displayed
    %----------------------------------------------------------------------
    
    % Voronoi cell colors/transparency
    fSegColors = 0.8 .* ones( size(c,1), 3 );
    fSegAlpha = zeros( size(c,1), 1 );
    
    % A list of nodes that will be displayed
    dispNodes = false(size(c,1), 1);
    
    for k = 1:numLabels
        
        % The nodes with the current label at the correct time point
        [ ~, locNode ] = ismember( labelIDx{k}, nodeIDx );
        locNode( locNode == 0 ) = [];
        
        if ~isempty( locNode )
            
            fSegColors( locNode, : ) = ...
                repmat( segColors(k, :), numel(locNode), 1 );
            
            if ismember( segLabels{k}, badLabels )
                fSegAlpha( locNode ) = 1;
            else
                fSegAlpha( locNode ) = 1;
                dispNodes( locNode ) = true;
            end

        end

    end
    
    %----------------------------------------------------------------------
    % Process Nested Patch Vertices/Colors for Generation Display
    %----------------------------------------------------------------------
    
    % Convert list of displayed nodes from logical vector to index list
    dispNodes = find(dispNodes);
    
    % The face connectivity list of the nested patch
    genFace = nan(numel(dispNodes), maxFaceSize);
    
    % The vertex coordinate list of the nested patch
    genVertex = zeros( sum(~isnan(voronoiFace(dispNodes, :)), 'all'), 2 );
    
    % The scale factor for the nested patch vertices
    nestedScale = 0.5;
    
    % Calculate the vertex positions
    vertexCount = 1;
    for k = 1:size(genFace, 1)
        
        curC = c{dispNodes(k)};
        
        curV = v(curC, :);
        curV = nestedScale .* (curV - mean(curV,1)) + mean(curV,1);
        
        genFace(k, 1:numel(curC)) = ...
            vertexCount:(vertexCount+numel(curC)-1);
        
        genVertex(vertexCount:(vertexCount+numel(curC)-1), :) = curV;
        
        vertexCount = vertexCount + numel(curC );
        
    end
    
    % Nested Voronoi cell colors
    fGenColors = genColors(curGen(dispNodes), :);
        
    %----------------------------------------------------------------------
    % Generate Images
    %----------------------------------------------------------------------
    fig = figure('units', 'normalized', 'outerposition', [0 0 1 1], ...
        'visible', 'off');
    
    % Basic background image
    hIm = imshow(cat(3, MIP, MIP, MIP));
    fig.OuterPosition = [0 0 1 1];
    
    hold on
    
    % The primary Voronoi tesselation with segement label colors
    patch( 'Faces', voronoiFace, 'Vertices', v, ...
        'FaceVertexCData', fSegColors, 'FaceColor', 'flat', ...
        'FaceVertexAlphaData', fSegAlpha, 'FaceAlpha', 'flat', ...
        'AlphaDataMapping', 'none', 'EdgeColor', 'k', 'LineWidth', 3 );
    
    % The nested patch showing cell cycle generation
    patch( 'Faces', genFace, 'Vertices', genVertex, ...
        'FaceVertexCData', fGenColors, 'FaceColor', 'flat', ...
        'FaceVertexAlphaData', 1, 'FaceAlpha', 'flat', ...
        'AlphaDataMapping', 'none', 'EdgeColor', 'none' );
    
    hold off
    
    set( gca, 'YDir', 'normal' );
    
    xlim(get(hIm, 'XData'));
    ylim(get(hIm, 'YData'));
    xticks([]); yticks([]);
    
    % Write to file -------------------------------------------------------
    saveFileName = fullfile( saveDir, ...
        sprintf( 'Parasegments_T%04d.tif', t ) );
    
    % cdata = export_fig( saveFileName, '-tif', '-transparent' );
    % cdata = export_fig;
    
    F = getframe;
    cdata = F.cdata;
    cdata = imresize(cdata, [1000 1000]);
    imwrite( cdata, saveFileName, 'tif' );
    
    close(fig);
    
end
    
    
    
    






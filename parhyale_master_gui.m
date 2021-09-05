function [G, divStruct] = parhyale_master_gui(Gin, ImFileName)
%PARHYALE_MASTER_GUI An interactive GUI that allows the user to curate and
%manually correct cell tracks and lineages constructed via automatic
%methods.  Additional modes exist to generate and curate parasegment labels
%that are propagated through the cell tracking. 
%
%   INPUT PARAMETERS:
%
%       - Gin:      A MATLAB digraph object. The nodes of the digraph
%                   represent all of the cells from every time point.  The
%                   connectivity of the digraph contains information about
%                   cell tracks from time point to time point.
%
%       - ImFileName:   The name of the image file over which to
%                       display the tracked cells at each given
%                       time point.  It is assumed that each time
%                       point is saved separately with identical
%                       names up to an integer which delineated the
%                       time point.
%
%   OUTPUT PARAMETERS:
%
%       - G:        The MATLAB digraph object representing the manually
%                   curated tracks and lineages
%
%       - divStruct     A structure containing information about the
%                       division events derived from the tracking graph
%
%   by Dillon Cislo 06/03/2020

%==========================================================================
% INPUT PROCESSING
%==========================================================================

% The initial output graph is just a copy of the input graph
G = Gin;

% Check that the input graph has the necessary fields
reqVars = {'T', 'UPix', 'Segment', 'Generation' };
assert( all(ismember(reqVars, G.Nodes.Properties.VariableNames)), ...
    'Input digraph is missing required fields!' );

% Extract the time points from the graph
timePoints = unique(G.Nodes.T);

% Validate generation counts
if any(G.Nodes.Generation < 1)
    error('Invalid generation counts!');
end

% Validate file name convention
validateattributes(ImFileName, {'char'}, {'scalartext'});

% Generate division struct from tracking graph
divStruct = assemble_div_struct(Gin);

%==========================================================================
% GENERATE FIGURE WINDOW
%==========================================================================

% The length scale of each marker face in units based on the x-axis 
% of each displayed image.  Used to construct scatter plots in tracking
% mode
markerSize = 20;

% Create the figure window ------------------------------------------------
fig = figure('units', 'normalized', 'outerposition', [0 0 1 1], ...
    'HitTest', 'off', 'KeyPressFcn', @kp_Fcn);

maximize(fig);

% The figure window position
fPos = fig.Position;

% Create the subplot axes for the three plots -----------------------------
axesList = gobjects(3, 1);

% Tracking Mode: The axis containing the image of the parent time point
% Segment Mode: The axis containing the image of the previous time point
prevAxPos = [ fPos(1)+fPos(3)/4, fPos(2)+fPos(4)/2, ...
    fPos(3)/4, fPos(4)/2 ];
axesList(1) = subplot('Position', prevAxPos);

% Tracking Mode: The axis containing the image of the fused time points
% Segment Mode: The axis containing the image of the current time point
curAxPos = [ fPos(1)+fPos(3)/2, fPos(2), ...
    fPos(3)/2, fPos(4) ];
axesList(2) = subplot('Position', curAxPos);

% Tracking Mode: The axis containing the image of the child time point
% Segment Mode: The axis containing the image of the next time point
nextAxPos = [ fPos(1)+fPos(3)/4, fPos(2), ...
    fPos(3)/4, fPos(4)/2 ];
axesList(3) = subplot('Position', nextAxPos);

% Slave the axis limits of all plots
linkaxes( axesList );

T0 = timePoints(1:3);

% Create a handle to a zoom utility ---------------------------------------
hZoom = zoom;
set(hZoom, 'ActionPostCallback', @zoomCallBack);

%==========================================================================
% CREATE UI CONTROLS
%==========================================================================

%--------------------------------------------------------------------------
% Text Box Controls
%--------------------------------------------------------------------------

% A label for the time point slider
uicontrol( 'Parent', fig, 'Style', 'text', ...
    'Units', 'normalized', 'Position', [ 0 0.95 0.0625 0.025 ], ...
    'String', 'Time Point', 'BackgroundColor', fig.Color, ...
    'FontWeight', 'bold', 'FontSize', 13 );

% A label for the GUI mode pop-up menu
uicontrol( 'Parent', fig, 'Style', 'text', ...
    'Units', 'normalized', 'Position', [ 0 0.915 0.1 0.025 ], ...
    'String', 'Choose GUI Mode:', 'BackgroundColor', fig.Color, ...
    'FontWeight', 'bold', 'FontSize', 13 );

% A label for the tracking mode buttons
uicontrol( 'Parent', fig, 'Style', 'text', ...
    'Units', 'normalized', 'Position', [0 0.88 0.15 0.025 ], ...
    'String', '    Tracking Mode Functions', ...
    'BackgroundColor', fig.Color, 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold', 'FontSize', 13 );

% A label for the segment label mode buttons
uicontrol( 'Parent', fig, 'Style', 'text', ...
    'Units', 'normalized', 'Position', [0 0.64 0.185 0.025 ], ...
    'String', '    Segment Label Mode Functions', ...
    'BackgroundColor', fig.Color, 'FontWeight', 'bold', ...
    'FontSize', 13, 'HorizontalAlignment', 'left' );

% A label for the general function buttons
uicontrol( 'Parent', fig, 'Style', 'text', ...
    'Units', 'normalized', 'Position', [0 0.4515 0.15 0.025 ], ...
    'String', '    General Functions', ...
    'BackgroundColor', fig.Color, 'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold', 'FontSize', 13 );

% A label for the generation type listbox
uicontrol( 'Parent', fig, 'Style', 'text', ...
    'Units', 'normalized', 'Position', [1/300 0.32 0.12 0.025], ...
    'String', 'Generations', 'BackgroundColor', fig.Color, ...
    'FontWeight', 'bold', 'FontSize', 13 );

% A label for the segment type listbox
uicontrol( 'Parent', fig, 'Style', 'text', ...
    'Units', 'normalized', 'Position', [(.12+2/300) 0.32 0.12 0.025], ...
    'String', 'Segment Labels', 'BackgroundColor', fig.Color, ...
    'FontWeight', 'bold', 'FontSize', 13 );

% A label for the time point of the first axis
axLabel1 = uicontrol( 'Parent', fig, 'Style', 'text', ...
    'Units', 'normalized', 'Position', [ 0.255 0.95 0.05 0.025 ], ...
    'String', sprintf('T = %d', T0(1)), 'BackgroundColor', 'k', ...
    'ForegroundColor', 'w', 'FontWeight', 'bold', 'FontSize', 13 );

% A label for the time point of the second axis
axLabel2 = uicontrol( 'Parent', fig, 'Style', 'text', ...
    'Units', 'normalized', 'Position', [ 0.505 0.95 0.05 0.025 ], ...
    'String', sprintf('T = %d', T0(2)), 'BackgroundColor', 'k', ...
    'ForegroundColor', 'w', 'FontWeight', 'bold', 'FontSize', 13 );

% A label for the time point of the third axis
axLabel3 = uicontrol( 'Parent', fig, 'Style', 'text', ...
    'Units', 'normalized', 'Position', [ 0.255 0.4625 0.05 0.025 ], ...
    'String', sprintf('T = %d', T0(3)), 'BackgroundColor', 'k', ...
    'ForegroundColor', 'w', 'FontWeight', 'bold', 'FontSize', 13 );

% A prompt for communicating with the user
userPrompt = uicontrol( 'Parent', fig, 'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [ 1/300 1/300 (0.25-2/300) (0.15-2/300) ], ...
    'String', [ 'Welcome to the Parhyale Cell Tracking/' ...
    'Parasegment Label Master GUI! ' ...
    'Navigate time points using the left/right arrow keys. ' ...
    'Switch between modes using the up/down arrow keys' ], ...
    'BackgroundColor', fig.Color, 'FontWeight', 'bold', 'FontSize', 12);

%--------------------------------------------------------------------------
% Button Controls
%--------------------------------------------------------------------------

% Tracking Mode Buttons ---------------------------------------------------

% Update parent node button
updateParentBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Update Parent', 'FontWeight', 'bold', 'FontSize', 12, ...
    'Units', 'normalized', ...
    'Position', [ 1/300 0.8345 0.12 0.045 ], ...
    'CallBack', @updateParent );

% Make progenitor button
makeProgenitorBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Make Progenitor', 'FontWeight', 'bold', 'FontSize', 12, ...
    'Units', 'normalized', ...
    'Position', [ (.12+2/300) 0.8345 0.12 0.045 ], ...
    'CallBack', @makeProgenitor );

% Add node button
addNodeBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Add Child Node', 'FontWeight', 'bold', 'FontSize', 12, ...
    'Units', 'normalized', ...
    'Position', [ 1/300 0.783 0.12 0.045 ], ...
    'CallBack', @addNode );

% Merge child nodes button
mergeChildBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Merge Child Nodes', 'FontWeight', 'bold', ...
    'FontSize', 12, 'Units', 'normalized', ...
    'Position', [ (.12+2/300) 0.783 0.12 0.045 ], ...
    'CallBack', @mergeChild );

% Remove parent button
RemoveParentBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Remove Parent Node', 'FontWeight', 'bold', ...
    'FontSize', 12, 'Units', 'normalized', ...
    'Position', [ 1/300 0.7315 0.12 0.045 ], ...
    'CallBack', @removeParent );

% Remove child button
removeChildBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Remove Child Node', 'FontWeight', 'bold', ...
    'FontSize', 12, 'Units', 'normalized', ...
    'Position', [ (.12+2/300) 0.7315 0.12 0.045 ], ...
    'CallBack', @removeChild );

% Update point matching button
updatePointMatchBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Update Point Matching', 'FontWeight', 'bold', ...
    'FontSize', 12, 'Units', 'normalized', ...
    'Position', [ 1/300 0.68 0.12 0.045 ], ...
    'CallBack', @updatePointMatch );

% Perform Demons-augmented point matching on the current time point
demonsBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Demons Matching', ...
    'FontWeight', 'bold', 'FontSize', 12', ...
    'Units', 'normalize', ...
    'Position', [ (.12+2/300) 0.68 0.12 0.045 ], ...
    'CallBack', @demonsTracking );

% Segment Label Mode Buttons ----------------------------------------------

% Set label button
setLabelBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Set Cell Label', ...
    'FontWeight', 'bold', 'FontSize', 12, ...
    'Units', 'normalized', ...
    'Position', [ 1/300 0.5945 0.12 0.045 ], ...
    'CallBack', @setLabel );

% Unset label button
unsetLabelBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Unset Cell Label', ...
    'FontWeight', 'bold', 'FontSize', 12, ...
    'Units', 'normalized', ...
    'Position', [ (.12+2/300) 0.5945 0.12 0.045 ], ...
    'CallBack', @unsetLabel );

% Add new label button
addLabelBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Add Label', ...
    'FontWeight', 'bold', 'FontSize', 12, ...
    'Units', 'normalized', ...
    'Position', [ 1/300 0.543 0.12 0.045 ], ...
    'CallBack', @addLabel );

% Remove label button
removeLabelBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Remove Label', ...
    'FontWeight', 'bold', ...
    'FontSize', 12, 'Units', 'normalized', ...
    'Position', [ (.12+2/300) 0.543 0.12 0.045 ], ...
    'CallBack', @removeLabel );

% Set generation button
setGenBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Set Generation', ...
    'FontWeight', 'bold', 'FontSize', 12, ...
    'Units', 'normalized', ...
    'Position', [ 1/300 0.4915 0.12 0.045 ], ...
    'CallBack', @setGen );

% Propagate generation  button
propGenBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Propagate Generations', ...
    'FontWeight', 'bold', 'FontSize', 12, ...
    'Units', 'normalized', ...
    'Position', [ (.12+2/300) 0.4915 0.12 0.045 ], ...
    'CallBack', @propGen );

% General Function Buttons ------------------------------------------------

% Reset zoom button
resetZoomBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Reset Zoom', ...
    'FontWeight', 'bold', ...
    'FontSize', 12, 'Units', 'normalized', ...
    'Position', [ 1/300 0.406 0.12 0.045 ], ...
    'CallBack', @resetZoom );

% Exit guide button
exitBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Exit GUI', ...
    'FontWeight', 'bold', ...
    'FontSize', 12, 'Units', 'normalized', ...
    'Position', [ (.12+2/300) 0.406 0.12 0.045 ], ...
    'CallBack', @exitGuide );

% Save Progress Btn
saveBtn = uicontrol( 'Parent', fig, 'Style', 'pushbutton', ...
    'String', 'Save Progress', 'FontWeight', 'bold', ...
    'FontSize', 12, 'Units', 'normalized', ...
    'Position', [ 1/300 0.3545 0.12 0.045 ], ...
    'CallBack', @saveProgress );

%--------------------------------------------------------------------------
% Toggle Button Controls
%--------------------------------------------------------------------------

% A general purpose user action button
userActionBtn = uicontrol( 'Parent', fig, 'Style', 'togglebutton', ...
    'Min', 0, 'Max', 1, 'Value', 0, ...
    'FontWeight', 'bold', 'FontSize', 12', ...
    'units', 'normalize', ...
    'Position', [ (.12+2/300) 0.3545 0.12 0.045 ] );

%--------------------------------------------------------------------------
% Slider Controls
%--------------------------------------------------------------------------

% A slider to help choose the current (child) time point
timeSld = uicontrol( 'Parent', fig, 'Style', 'slider', 'Value', 2, ...
    'Units', 'normalized', 'Position', [ 0.0625 0.95 0.125 0.025 ], ...
    'Min', 1, 'Max', numel(timePoints), ...
    'SliderStep', [1 1]/(numel(timePoints)-1), ...
    'CallBack', @updateTimeSld );

%--------------------------------------------------------------------------
% Edit Controls
%--------------------------------------------------------------------------

% An edit to help choose the current (child) time point
timeEdit = uicontrol( 'Parent', fig, 'Style', 'edit', ...
    'Units', 'normalized', ...
    'String', num2str(timePoints(timeSld.Value)), ..., 
    'Position', [ (0.1875 + 0.00625) 0.95 (0.8 * 0.0625) 0.025 ], ...
    'BackgroundColor', fig.Color', 'CallBack', @updateTimeEdit );

%--------------------------------------------------------------------------
% Listbox Controls
%--------------------------------------------------------------------------
genBox = uicontrol( 'Parent', fig, 'Style', 'listbox', ...
    'Units', 'normalized', 'FontSize', 15, 'FontWeight', 'bold', ...
    'Position', [1/300 (0.15+1/300) (0.125-2/300) (0.17-2/300) ], ...
    'String', [], 'BackgroundColor', 'w' );

labelBox = uicontrol( 'Parent', fig, 'Style', 'listbox', ...
    'Units', 'normalized', 'FontSize', 15, 'FontWeight', 'bold', ...
    'Position', ...
    [ (0.125+1/300) (0.15+1/300) (0.125-2/300) (0.17-2/300) ], ...
    'String', [], 'BackgroundColor', 'w' );

% Some hacky list box layout features -------------------------------------

jScrollPane_Gen = java(findjobj(genBox));
jListBox_Gen = jScrollPane_Gen.getViewport.getView;
set( jListBox_Gen, 'LayoutOrientation', 1, 'VisibleRowCount', 21 );

jScrollPane_Label = java(findjobj(labelBox));
jListbox_Label = jScrollPane_Label.getViewport.getView;
set( jListbox_Label, 'LayoutOrientation', 1, 'VisibleRowCount', 21 );

% Set the initial label box contents --------------------------------------

% The unique parasegment labels
segLabels = unique( G.Nodes.Segment );
segLabels( strcmp(segLabels, 'none') ) = [];

% Determine which nodes should contribute to the generation coloring scheme
has_Label = ismember(G.Nodes.Segment, segLabels);

% The number of unique generations
num_Gen = max(G.Nodes(has_Label,:).Generation);

% Assign a unique color to each label and generation
seg_Colors = 255 .* distinguishable_colors( numel(segLabels) + num_Gen, ...
    [0 0 0; 1 1 1] );

gen_Colors = seg_Colors(1:num_Gen, :);
seg_Colors = seg_Colors((num_Gen+1):end, :);

% Set the generation string
gen_Types = { 'PSPR', 'Post-Wave 1', ...
    'Post-Wave 2', 'Diff Cleavage %d'};
gen_String = cell(size(gen_Colors,1),1);
for ii = 1:size(gen_Colors,1)
    
    if (ii < 4)
        gen_String{ii} = sprintf( ...
            ['<HTML><BODY bgcolor = "rgb(%f,%f,%f">' ...
            gen_Types{ii} '</BODY></HTML>'], ...
            gen_Colors(ii,1), gen_Colors(ii,2), gen_Colors(ii,3) );
    else
        gen_String{ii} = sprintf( ...
            ['<HTML><BODY bgcolor = "rgb(%f,%f,%f">' ...
            gen_Types{4} '</BODY></HTML>'], ...
            gen_Colors(ii,1), gen_Colors(ii,2), gen_Colors(ii,3), ii-3 );
    end
    
end

genBox.String = gen_String;
        
% Set the label string
label_String = cell(size(segLabels));
for ii = 1:numel(segLabels)
    label_String{ii} = sprintf( ...
        ['<HTML><BODY bgcolor="rgb(%f,%f,%f)">' ...
        segLabels{ii} '</BODY></HTML>'], ...
        seg_Colors(ii,1), seg_Colors(ii,2), seg_Colors(ii,3) );
end

labelBox.String = label_String;

%--------------------------------------------------------------------------
% Pop-Up Menu Controls
%--------------------------------------------------------------------------

% A pop-up menu to select the mode in which the GUI operates
modeMenu = uicontrol( 'Parent', fig, 'Style', 'popupmenu', ...
    'Units', 'normalized', ...
    'String', {'Tracking Mode', 'Segment Label Mode'}, ...
    'Position', [ 0.105 0.915 0.12 0.025 ], ...
    'FontWeight', 'bold', 'FontSize', 13, ...
    'Value', 2, 'Callback', @setGUIMode );

%==========================================================================
% DRAW FIGURES
%==========================================================================
draw_seg( G, T0, segLabels, axesList, ImFileName);

%==========================================================================
uiwait(fig); % Prevent Function Output Until Correction is Completed
%==========================================================================

%==========================================================================
%**************************************************************************
%==========================================================================
% CALLBACK FUNCTIONS
%==========================================================================
%**************************************************************************
%==========================================================================

%==========================================================================
% TRACKING MODE FUNCTIONS
%==========================================================================

%--------------------------------------------------------------------------
% UPDATE PARENT
%--------------------------------------------------------------------------
% This function allows replaces the parent node of a user selected child
% node with another user selected node from the previous time point

    function updateParent( ~, ~, ~ )
        
        if (modeMenu.Value == 2)
            
            userPrompt.String = ['The "Update Parent" function ' ...
                'cannot be used in segment label mode'];
            
            return;
            
        end
        
        % Extract data about the current/previous time points -------------
        
        cur_T = str2double(timeEdit.String); % The current (child) time
        tidx = find(timePoints == cur_T); % The ID of that time point
        prev_T = timePoints(tidx-1); % The previous (parent) time
        
        % The node IDs of nodes at the current (child) time
        curIDx = find( G.Nodes.T == cur_T );
        
        % The node IDs of nodes at the previous (parent) time
        prevIDx = find( G.Nodes.T == prev_T );
        
        % The (x,y)-coordinates of all child nodes in pixel space
        curU = G.Nodes(curIDx,:).UPix;
        
        % The (x,y)-coordinates of all parent nodes
        prevU = G.Nodes(prevIDx,:).UPix;
        
        XLim = get(gca, 'XLim'); % The x-axis limits
        YLim = get(gca, 'YLim'); % The y-axis limits
        
        % Select the child node -------------------------------------------
        
        % Prompt for user input
        userPrompt.String = 'Please select one child node';
        [xC, yC] = ginput(1);
        
        % Check the the specified point lies within the axis limits
        if ( (xC <= XLim(1)) || (XLim(2) <= xC) || ...
                (yC <= YLim(1)) || (YLim(2) <= yC) )
            
            userPrompt.String = [ 'Please choose a point that '...
                'lies within the axes limits!' ];
            
            return;
            
        end
        
        % Point match to find the node ID of the child
        childID = curU - repmat([xC, yC], size(curU,1), 1);
        childID = sqrt( sum( childID.^2, 2 ) );
        [~, childID] = min(childID);
        childID = curIDx(childID);
        
        % Display the chosen child point
        subplot(axesList(2));
        hold on
        hCirc = viscircles( G.Nodes(childID,:).UPix, ...
            markerSize/1.5, 'Color', 'w', 'LineWidth', markerSize/8 );
        hold off
        
        % Select the new parent node --------------------------------------
        
        % Prompt for user input
        userPrompt.String = 'Please select one new parent node';
        [xP, yP] = ginput(1);
        
        % Check the the specified point lies within the axis limits
        if ( (xP <= XLim(1)) || (XLim(2) <= xP) || ...
                (yP <= YLim(1)) || (YLim(2) <= yP) )
            
            userPrompt.String = [ 'Please choose a point that '...
                'lies within the axes limits!' ];
            
            delete(hCirc);
            
            return;
            
        end
        
        % Point match to find the node ID of the new parent
        parentID = prevU - repmat([xP, yP], size(prevU,1), 1);
        parentID = sqrt( sum( parentID.^2, 2 ) );
        [~, parentID] = min(parentID);
        parentID = prevIDx(parentID);
        
        % Update the edge in the digraph to reflect the new parent --------
        
        % The ID of input edge to the child (there should be only one)
        eID = find(G.Edges.EndNodes(:,2) == childID);
        
        if (numel(eID) > 1)
            
            userPrompt.String = ['Selected child node has multiple ' ...
                'parents'];
            
            return;
            
        end
        
        % Remove the old edge
        if ~isempty(eID),  G = rmedge(G, eID); end 
        
        G = addedge(G, parentID, childID); % Add the new edge
        
        % Update the division information ---------------------------------
        divStruct = assemble_div_struct(G);
        
        % Change the label properties -------------------------------------
        
        % The child node inherits the segment label of the parent
        G.Nodes(childID,:).Segment = G.Nodes(parentID,:).Segment;
        
        % The child node inherits the generation label of the parent
        % accounting for the possibility of a division event
        G.Nodes(childID,:).Generation = ...
            G.Nodes(parentID,:).Generation + ...
            ismember(parentID, divStruct.parentID);
        
        % Propagate label properties to descendants -----------------------
        
        % Extract all descendants of the current node
        descNodes = dfsearch( G, childID, ...
            {'discovernode', 'edgetonew' } );
        
        % Extract edges defining lineage of the current node
        descEdges = descNodes.Edge;
        descEdges( any(isnan(descEdges),2), : ) = [];
        
        descSources = descEdges(:,1);
        descSinks = descEdges(:,2);
        
        % Determine which edge source are division events
        isDiv = ismember(descSources, divStruct.parentID);
        
        % Update segment labels for descendants
        G.Nodes(descSinks,:).Segment = ...
            repmat(G.Nodes(childID,:).Segment, numel(descSinks), 1);
        
        % Update generation labels for descendants
        for j = 1:size(descEdges, 1)
            
            if isDiv(j)
                G.Nodes(descSinks(j),:).Generation = ...
                    G.Nodes(descSources(j),:).Generation + 1;
            else
                G.Nodes(descSinks(j),:).Generation = ...
                    G.Nodes(descSources(j),:).Generation;
            end
            
        end
        
        % Re-draw display images ------------------------------------------
        T = [prev_T, cur_T, -1];
        
        draw_tracks( G, T, axesList, ImFileName, markerSize);
        
        set(axesList(2), 'XLim', XLim);
        set(axesList(2), 'YLim', YLim);
        
        scaleFigChildren(fig, markerSize);
        
    end

%--------------------------------------------------------------------------
% MAKE PROGENITOR
%--------------------------------------------------------------------------
% This function makes a given node at the child time point a progenitor,
% i.e. removes all edges with the given node as the target

    function makeProgenitor( ~, ~, ~ )
        
        % Extract data about the current/previous time points -------------
        
        cur_T = str2double(timeEdit.String); % The current (child) time
        tidx = find(timePoints == cur_T); % The ID of that time point
        prev_T = timePoints(tidx-1); % The previous (parent) time
        
        % The node IDs of nodes at the current (child) time
        curIDx = find( G.Nodes.T == cur_T );
        
        % The (x,y)-coordinates of all child nodes in pixel space
        curU = G.Nodes(curIDx,:).UPix;
        
        XLim = get(gca, 'XLim'); % The x-axis limits
        YLim = get(gca, 'YLim'); % The y-axis limits
        
        % Select the child node -------------------------------------------
        
        % Prompt for user input
        userPrompt.String = 'Please select one child node';
        [xC, yC] = ginput(1);
        
        % Check the the specified point lies within the axis limits
        if ( (xC <= XLim(1)) || (XLim(2) <= xC) || ...
                (yC <= YLim(1)) || (YLim(2) <= yC) )
            
            userPrompt.String = [ 'Please choose a point that '...
                'lies within the axes limits!' ];
            
            return;
            
        end
        
        % Point match to find the node ID of the child
        childID = curU - repmat([xC, yC], size(curU,1), 1);
        childID = sqrt( sum( childID.^2, 2 ) );
        [~, childID] = min(childID);
        childID = curIDx(childID);
        
        % Find and remove all edges with the current node as a target -----
        
        % The ID of input edge to the child (there should be only one)
        eID = find(G.Edges.EndNodes(:,2) == childID);
        
        G = rmedge(G, eID); % Remove the old edge
        
        % Update division information -------------------------------------
        divStruct = assemble_div_struct(G);
        
        % Change the label properties -------------------------------------
        % Note the segment label of the progenitor node remains unchanged
        % after being made a progenitor node
        
        % Set the generation label of the progenitor node to one
        G.Nodes(childID,:).Generation = 1;
        
        % Propagate label properties to descendants -----------------------
        
        % Extract all descendants of the current node
        descNodes = dfsearch( G, childID, ...
            {'discovernode', 'edgetonew' } );
        
        % Extract edges defining lineage of the current node
        descEdges = descNodes.Edge;
        descEdges( any(isnan(descEdges),2), : ) = [];
        
        descSources = descEdges(:,1);
        descSinks = descEdges(:,2);
        
        % Determine which edge source are division events
        isDiv = ismember(descSources, divStruct.parentID);
        
        % Update segment labels for descendants
        G.Nodes(descSinks,:).Segment = ...
            repmat(G.Nodes(childID,:).Segment, numel(descSinks), 1);
        
        % Update generation labels for descendants
        for j = 1:size(descEdges, 1)
            
            if isDiv(j)
                G.Nodes(descSinks(j),:).Generation = ...
                    G.Nodes(descSources(j),:).Generation + 1;
            else
                G.Nodes(descSinks(j),:).Generation = ...
                    G.Nodes(descSources(j),:).Generation;
            end
            
        end
        
        % Re-draw display images ------------------------------------------
        T = [prev_T, cur_T, -1];
        if (tidx < numel(timePoints)), T(3) = timePoints(tidx+1); end
        
        if (modeMenu.Value == 1)
            
            draw_tracks( G, T, axesList, ImFileName, markerSize);
            
            set(axesList(2), 'XLim', XLim);
            set(axesList(2), 'YLim', YLim);
            
            scaleFigChildren(fig, markerSize);
            
        else
            
            draw_seg( G, T, segLabels, axesList, ImFileName );
            
            set(axesList(2), 'XLim', XLim);
            set(axesList(2), 'YLim', YLim);
            
        end
        
    end

%--------------------------------------------------------------------------
% ADD NODE
%--------------------------------------------------------------------------
% This function allows the user to add a child node at the current time at
% a specified location with a given parent node

    function addNode( ~, ~, ~ )
        
        if (modeMenu.Value == 2)
            
            userPrompt.String = ['The "Add Node" function cannot ' ...
                'be used in segment label mode'];
            
            return;
            
        end
        
        % Extract data about the current/previous time points -------------
        
        cur_T = str2double(timeEdit.String); % The current (child) time
        tidx = find(timePoints == cur_T); % The ID of that time point
        prev_T = timePoints(tidx-1); % The previous (parent) time
        
        % The node IDs of nodes at the previous (parent) time
        prevIDx = find( G.Nodes.T == prev_T );
        
        % The (x,y)-coordinates of all parent nodes
        prevU = G.Nodes(prevIDx,:).UPix;
        
        XLim = get(gca, 'XLim'); % The x-axis limits
        YLim = get(gca, 'YLim'); % The y-axis limits
        
        % Select the location for the new child node ----------------------
        
        % Prompt for user input
        userPrompt.String = [ 'Please choose a location for ' ...
            'the new child node' ];
        [xC, yC] = ginput(1);
        
        % Check the the specified point lies within the axis limits
        if ( (xC <= XLim(1)) || (XLim(2) <= xC) || ...
                (yC <= YLim(1)) || (YLim(2) <= yC) )
            
            userPrompt.String = [ 'Please choose a point that '...
                'lies within the axes limits!' ];
            
            return;
            
        end
        
        
        % Select the new parent node --------------------------------------
        
        % Prompt for user input
        userPrompt.String = 'Please select one new parent node';
        [xP, yP] = ginput(1);
        
        % Check the the specified point lies within the axis limits
        if ( (xP <= XLim(1)) || (XLim(2) <= xP) || ...
                (yP <= YLim(1)) || (YLim(2) <= yP) )
            
            userPrompt.String = [ 'Please choose a point that '...
                'lies within the axes limits!' ];
            
            return;
            
        end
        
        % Point match to find the node ID of the new parent
        parentID = prevU - repmat([xP, yP], size(prevU,1), 1);
        parentID = sqrt( sum( parentID.^2, 2 ) );
        [~, parentID] = min(parentID);
        parentID = prevIDx(parentID);
        
        % Update the graph structure --------------------------------------
        
        % The new node will inherit the label of its parent
        newLabel = G.Nodes(parentID,:).Segment;
        
        % Determine the cell cycle generation of the new node
        newGen = G.Nodes(parentID,:).Generation + ...
            isempty(successors(G, parentID));
        
        % Create a node table for the new node
        nodeTable = G.Nodes(1,:);
        
        % Extract the names of the node variables
        varNames = G.Nodes.Properties.VariableNames;
        
        % Nullify elements in the table
        % NOTE: These are NOT all possible data types that could be
        % contained in the table
        for i = 1:size(nodeTable,2)
            
            if iscell(nodeTable.(varNames{i}))
                
                if ischar(nodeTable.(varNames{i}){1})
                    
                    nodeTable.(varNames{i}) = ...
                        repmat({''}, size(nodeTable.(varNames{i})));
                    
                else
                    
                    nodeTable.(varNames{i}) = ...
                        cell(size(nodeTable.(varNames{i})));
                    
                end
            
            elseif isnumeric(nodeTable.(varNames{i}))
                
                nodeTable.(varNames{i}) = ...
                    zeros(size(nodeTable.(varNames{i})), ...
                    class(nodeTable.(varNames{i})));
                
            elseif isnan(nodeTable.(varNames{i}))
                
                nodeTable.(varNames{i}) = ...
                    nan(size(nodeTable.(varNames{i})), ...
                    class(nodeTable.(varNames{i})));
                
            elseif islogical(nodeTable.(varNames{i}))
                
                nodeTable.(varNames{i}) = ...
                    false(size(nodeTable.(varNames{i})));
                
            end
            
        end

        nodeTable.UPix = [xC, yC];
        nodeTable.T = cur_T;
        nodeTable.Segment = newLabel;
        nodeTable.Generation = newGen;
        
        % Add the new node to the graph
        G = addnode( G, nodeTable );
        
        % Add the incoming edge for the new node
        G = addedge( G, parentID, numnodes(G) );
        
        % Update division information -------------------------------------
        divStruct = assemble_div_struct(G);
        
        % Re-draw display images ------------------------------------------
        T = [prev_T, cur_T, -1];
        
        draw_tracks( G, T, axesList, ImFileName, markerSize);
        
        set(axesList(2), 'XLim', XLim);
        set(axesList(2), 'YLim', YLim);
        
        scaleFigChildren(fig, markerSize);
        
    end

%--------------------------------------------------------------------------
% REMOVE PARENT
%--------------------------------------------------------------------------
% This function allows the user to remove a parent node from the tracking
% structure

    function removeParent( ~, ~, ~ )
        
        if (modeMenu.Value == 2)
            
            userPrompt.String = ['The "Remove Parent" function ' ...
                'cannot be used in segment label mode'];
            
            return;
            
        end
        
        % Extract data about the current/previous time points -------------
        
        cur_T = str2double(timeEdit.String); % The current (child) time
        tidx = find(timePoints == cur_T); % The ID of that time point
        prev_T = timePoints(tidx-1); % The previous (parent) time
        
        % The node IDs of nodes at the previous (parent) time
        prevIDx = find( G.Nodes.T == prev_T );
        
        % The (x,y)-coordinates of all parent nodes
        prevU = G.Nodes(prevIDx,:).UPix;
        
        XLim = get(gca, 'XLim'); % The x-axis limits
        YLim = get(gca, 'YLim'); % The y-axis limits
        
        % Select the new parent node --------------------------------------
        
        % Prompt for user input
        userPrompt.String = 'Please select a parent node to remove';
        [xP, yP] = ginput(1);
        
        % Check the the specified point lies within the axis limits
        if ( (xP <= XLim(1)) || (XLim(2) <= xP) || ...
                (yP <= YLim(1)) || (YLim(2) <= yP) )
            
            userPrompt.String = [ 'Please choose a point that '...
                'lies within the axes limits!' ];
            
            return;
            
        end
        
        % Point match to find the node ID of the new parent
        parentID = prevU - repmat([xP, yP], size(prevU,1), 1);
        parentID = sqrt( sum( parentID.^2, 2 ) );
        [~, parentID] = min(parentID);
        parentID = prevIDx(parentID);
        
        % Update the labels of any children and their descendants ---------
        % The children will retain their segment labels
        
        childIDx = successors(G, parentID);
        if ~isempty(childIDx)
            
            % Make all children progenitors
            G.Nodes(childIDx,:).Generation = ones(numel(childIDx), 1);
            
            % Propagate label properties to descendants -------------------
            
            for i = 1:numel(childIDx)
                
                childID = childIDx(i);
                
                % Extract all descendants of the current node
                descNodes = dfsearch( G, childID, ...
                    {'discovernode', 'edgetonew' } );
                
                % Extract edges defining lineage of the current node
                descEdges = descNodes.Edge;
                descEdges( any(isnan(descEdges),2), : ) = [];
                
                descSources = descEdges(:,1);
                descSinks = descEdges(:,2);
                
                % Determine which edge source are division events
                isDiv = ismember(descSources, divStruct.parentID);
                
                % Update segment labels for descendants
                G.Nodes(descSinks,:).Segment = ...
                    repmat(G.Nodes(childID,:).Segment, ...
                    numel(descSinks), 1);
                
                % Update generation labels for descendants
                for j = 1:size(descEdges, 1)
                    
                    if isDiv(j)
                        G.Nodes(descSinks(j),:).Generation = ...
                            G.Nodes(descSources(j),:).Generation + 1;
                    else
                        G.Nodes(descSinks(j),:).Generation = ...
                            G.Nodes(descSources(j),:).Generation;
                    end
                    
                end
                
            end
            
        end
        
        % Remove the node from the tracking structure ---------------------
        G = rmnode(G, parentID);
        
        % Update division information -------------------------------------
        divStruct = assemble_div_struct(G);
        
        % Re-draw display images ------------------------------------------
        T = [prev_T, cur_T, -1];
        
        draw_tracks( G, T, axesList, ImFileName, markerSize);
        
        set(axesList(2), 'XLim', XLim);
        set(axesList(2), 'YLim', YLim);
        
        scaleFigChildren(fig, markerSize);
        
    end

%--------------------------------------------------------------------------
% REMOVE CHILD
%--------------------------------------------------------------------------
% This function allows the user to remove a child node from the tracking
% structure

    function removeChild( ~, ~, ~ )
        
        % Extract data about the current/previous time points -------------
        
        cur_T = str2double(timeEdit.String); % The current (child) time
        tidx = find(timePoints == cur_T); % The ID of that time point
        prev_T = timePoints(tidx-1); % The previous (parent) time
        
        % The node IDs of nodes at the current (child) time
        curIDx = find( G.Nodes.T == cur_T );
        
        % The (x,y)-coordinates of all child nodes
        curU = G.Nodes(curIDx,:).UPix;
        
        XLim = get(gca, 'XLim'); % The x-axis limits
        YLim = get(gca, 'YLim'); % The y-axis limits
        
        % Select the child node -------------------------------------------
        
        % Prompt for user input
        userPrompt.String = 'Please select one child node';
        [xC, yC] = ginput(1);
        
        % Check the the specified point lies within the axis limits
        if ( (xC <= XLim(1)) || (XLim(2) <= xC) || ...
                (yC <= YLim(1)) || (YLim(2) <= yC) )
            
            userPrompt.String = [ 'Please choose a point that '...
                'lies within the axes limits!' ];
            
            return;
            
        end
        
        % Point match to find the node ID of the child
        childID = curU - repmat([xC, yC], size(curU,1), 1);
        childID = sqrt( sum( childID.^2, 2 ) );
        [~, childID] = min(childID);
        childID = curIDx(childID);
        
        % Update the label of any grandchildren and their descendants -----
        % The grandchildren will retain their segment labels
        
        grandChildIDx = successors(G, childID);
        if ~isempty(grandChildIDx)
            
            % Make all grandchildren progenitors
            G.Nodes(grandChildIDx,:).Generation = ...
                ones(numel(grandChildIDx), 1);
            
            % Propagate label properties to descendants -------------------
            
            for i = 1:numel(grandChildIDx)
                
                grandChildID = grandChildIDx(i);
                
                % Extract all descendants of the current node
                descNodes = dfsearch( G, grandChildID, ...
                    {'discovernode', 'edgetonew' } );
                
                % Extract edges defining lineage of the current node
                descEdges = descNodes.Edge;
                descEdges( any(isnan(descEdges),2), : ) = [];
                
                descSources = descEdges(:,1);
                descSinks = descEdges(:,2);
                
                % Determine which edge source are division events
                isDiv = ismember(descSources, divStruct.parentID);
                
                % Update segment labels for descendants
                G.Nodes(descSinks,:).Segment = ...
                    repmat(G.Nodes(grandChildID,:).Segment, ...
                    numel(descSinks), 1);
                
                % Update generation labels for descendants
                for j = 1:size(descEdges, 1)
                    
                    if isDiv(j)
                        G.Nodes(descSinks(j),:).Generation = ...
                            G.Nodes(descSources(j),:).Generation + 1;
                    else
                        G.Nodes(descSinks(j),:).Generation = ...
                            G.Nodes(descSources(j),:).Generation;
                    end
                    
                end
                
            end
            
        end
        
        % Remove node from tracking structure -----------------------------
        G = rmnode(G, childID);
        
        % Update division information -------------------------------------
        divStruct = assemble_div_struct(G);
        
        % Re-draw display images ------------------------------------------
        T = [prev_T, cur_T, -1];
        if (tidx < numel(timePoints)), T(3) = timePoints(tidx+1); end
        
        if (modeMenu.Value == 1)
            
            draw_tracks( G, T, axesList, ImFileName, markerSize);
            
            set(axesList(2), 'XLim', XLim);
            set(axesList(2), 'YLim', YLim);
            
            scaleFigChildren(fig, markerSize);
            
        else
            
            draw_seg( G, T, segLabels, axesList, ImFileName );
            
            set(axesList(2), 'XLim', XLim);
            set(axesList(2), 'YLim', YLim);
            
        end
        
    end

%--------------------------------------------------------------------------
% MERGE CHILD
%--------------------------------------------------------------------------
% This function allows the user to merge a set of nodes at the current 
% time.  All of the nodes being merged must share the same parent.

    function mergeChild( ~, ~, ~ )
        
        if (modeMenu.Value == 2)
            
            userPrompt.String = ['The "Merge Child" function ' ...
                'cannot be used in segment label mode'];
            
            return;
            
        end
        
        % Extract data about the current/previous/next time points --------
        
        cur_T = str2double(timeEdit.String); % The current (child) time
        tidx = find(timePoints == cur_T); % The ID of that time point
        prev_T = timePoints(tidx-1); % The previous (parent) time
        
        % The node IDs of nodes at the current (child) time
        curIDx = find( G.Nodes.T == cur_T );
        
        % The (x,y)-coordinates of all child nodes
        curU = G.Nodes(curIDx,:).UPix;
        
        % Determine if the current time point is the last time point
        isLastTime = tidx == numel(timePoints);
        
        if ~isLastTime
            % The ID of the next (grandchild) time point
            next_T = timePoints(tidx+1);
        else
            next_T = [];
        end
            
        XLim = get(gca, 'XLim'); % The x-axis limits
        YLim = get(gca, 'YLim'); % The y-axis limits
        
        % Modify user interaction strings 
        userPrompt.String = [ 'Please select child nodes to merge. ' ...
            'Push button to indicate when selection is complete' ];
        
        userActionBtn.String = 'Finish Selection';
        
        % Re-set the user action button in case the user forgot
        userActionBtn.Value = 0;
        
        % Extract child nodes to merge ------------------------------------
        
        % The node IDs to merge
        childIDx = []; 
        
        % A list of graphics handles to indicate the selected nodes
        hCircs = [];
        
        count = 0; % An indexing variable
        while true
            
            % The (x,y)-coordinates of the current selection
            [xC, yC] = ginput(1);
            drawnow()
            
            % Check the the specified point lies within the axis limits
            if ( (xC <= XLim(1)) || (XLim(2) <= xC) || ...
                    (yC <= YLim(1)) || (YLim(2) <= yC) )
 
                if userActionBtn.Value
                    break;
                else
                    userPrompt.String = [ 'Please choose a point that '...
                        'lies within the axes limits!' ];
                    continue;
                end

            end
            
            % Point match to find the node ID of the current selection
            childID = curU - repmat([xC, yC], size(curU,1), 1);
            childID = sqrt( sum( childID.^2, 2 ) );
            [~, childID] = min(childID);
            childID = curIDx(childID);
            
            if ismember(childID, childIDx)
                userPrompt.String = [ 'Selected node has already ' ...
                    'been picked. Please choose another' ];
                continue;
            end
            
            % Increment the index variable
            count = count + 1;
            
            % Add the current selection to the list of nodes to merge
            childIDx(count) = childID;
            
            % Visualize current selection on the fused axis
            subplot(axesList(2));
            hold on
            hCircs(count) = viscircles( G.Nodes(childID,:).UPix, ...
                markerSize/1.5, 'Color', 'w', 'LineWidth', markerSize/8 );
            hold off
            
        end
        
        % Re-set the user action button
        userActionBtn.Value = 0;
        
        % Delete selection displays for clarity
        for i = 1:numel(hCircs)
            delete(hCircs(i));
        end
        
        % Update interaction strings
        userActionBtn.String = [];
       
        if ~isempty(childIDx)
            userPrompt.String = 'Selection finished. Merging nodes ... ';
            drawnow
        else
            userPrompt.String = 'No nodes selected';
            return;
        end
        
        %------------------------------------------------------------------
        % Update tracking graph structure
        %------------------------------------------------------------------
        
        % Handle the shared parent ----------------------------------------
        
        % Find the ID of the shared parent node
        parentID = zeros(numel(childIDx), 1);
        for i = 1:numel(childIDx)
            parentID(i) = predecessors(G, childIDx(i));
        end
        
        parentID = unique(parentID);
        
        % Check that all selected nodes have the same parent
        if ( numel(parentID) ~= 1 )
            
            userPrompt.String = ['All selected nodes must share the ' ...
                'same parent in order to be merged.  Update parents ' ...
                'and try again' ];
            
            return;
            
        end
        
        % The (x,y)-coordinates of the parent node
        parentU = G.Nodes(parentID,:).UPix;
        
        % Handle the shared descendants of nodes being merged -------------
        
        grandChildIDx = []; % The node IDs of the descendants
        grandChildU = []; % The (x,y)-coordinates of the descendants
        
        if ~isLastTime
            
            for i = 1:numel(childIDx)
                gcIDx = successors(G, childIDx(i));
                grandChildIDx = [ grandChildIDx, ...
                    reshape( gcIDx, 1, numel(gcIDx) ) ];
            end
            
            if ~isempty(grandChildIDx)
                if ~isequal(grandChildIDx, unique(grandChildIDx, 'stable'))
                    
                    userPrompt.String = [ 'Descendants of selected ' ...
                        'nodes have more than one ancestor!' ];
                    
                    return;
                    
                end
            end
            
            grandChildU = G.Nodes(grandChildIDx,:).UPix;
            
        end
        
        % Find the coordinate of the new node -----------------------------
        % The coordinate of the new node will be the average of the
        % coordinate of the selected nodes
        UNew = mean(G.Nodes(childIDx, :).UPix, 1);
        
        % Node Handling ---------------------------------------------------
        
        % The new node will inherit the label of its parent
        newLabel = G.Nodes(parentID,:).Segment;
        
        % Determine the cell cycle generation of the new node
        newGen = G.Nodes(parentID,:).Generation + ...
            ismember(parentID, divStruct.parentID);
        
        % Create a node table for the new node
        nodeTable = G.Nodes(1,:);
        
        % Extract the names of the node variables
        varNames = G.Nodes.Properties.VariableNames;
        
        % Nullify elements in the table
        % NOTE: These are NOT all possible data types that could be
        % contained in the table
        for i = 1:size(nodeTable,2)
            
            if iscell(nodeTable.(varNames{i}))
                
                if ischar(nodeTable.(varNames{i}){1})
                    
                    nodeTable.(varNames{i}) = ...
                        repmat({''}, size(nodeTable.(varNames{i})));
                    
                else
                    
                    nodeTable.(varNames{i}) = ...
                        cell(size(nodeTable.(varNames{i})));
                    
                end
            
            elseif isnumeric(nodeTable.(varNames{i}))
                
                nodeTable.(varNames{i}) = ...
                    zeros(size(nodeTable.(varNames{i})), ...
                    class(nodeTable.(varNames{i})));
                
            elseif isnan(nodeTable.(varNames{i}))
                
                nodeTable.(varNames{i}) = ...
                    nan(size(nodeTable.(varNames{i})), ...
                    class(nodeTable.(varNames{i})));
                
            elseif islogical(nodeTable.(varNames{i}))
                
                nodeTable.(varNames{i}) = ...
                    false(size(nodeTable.(varNames{i})));
                
            end
            
        end

        nodeTable.UPix = UNew ;
        nodeTable.T = cur_T;
        nodeTable.Segment = newLabel;
        nodeTable.Generation = newGen;
        
        % Remove all of the old nodes
        G = rmnode(G, childIDx);
        
        % Add the new node to the graph
        G = addnode( G, nodeTable );
        
        % Edge handling ---------------------------------------------------
        % NOTE: the removal and addition of the merged nodes will have
        % changed the node IDs of the parent and descendants.  These IDs
        % must be updated before the edges can be properly added back into
        % the tracking structure
        
        % The IDs of all the nodes at the previous time
        prevIDx = find(G.Nodes.T == prev_T);
        
        % The (x,y)-coordinates of all the nodes at the previous time
        prevU = G.Nodes(prevIDx,:).UPix;
        
        % Point match to find the new ID of the parent node
        newParentID = prevU - repmat(parentU, numel(prevIDx), 1);
        newParentID = sqrt( sum( newParentID.^2, 2 ) );
        [~, newParentID] = min(newParentID);
        newParentID = prevIDx( newParentID );
        
        % Add the incoming edge
        G = addedge(G, newParentID, numnodes(G));
        
        if ~isLastTime && ~isempty(grandChildIDx)
            
            % The IDs of all the nodes at the next time
            nextIDx = find(G.Nodes.T == next_T);
            
            % The (x,y)-coordinates of all the nodes at the next time
            nextU = G.Nodes(nextIDx,:).UPix;
            
            % Point match to find the new ID of the descendants
            newGrandChildIDx = zeros(size(grandChildIDx));
            for i = 1:numel(grandChildIDx)
                
                gID = nextU - ...
                    repmat(grandChildU(i,:), numel(nextIDx), 1);
                gID = sqrt( sum( gID.^2, 2 ) );
                [~, gID] = min(gID);
                newGrandChildIDx(i) = nextIDx( gID );
                
            end
            
            % Add all outgoing edges
            G = addedge( G, ...
                repmat( numnodes(G), numel(grandChildIDx), 1 ), ...
                newGrandChildIDx );
            
        end
        
        % Update division information -------------------------------------
        divStruct = assemble_div_struct(G);
        
        % Propagate label properties to descendants -----------------------
        
        % Extract all descendants of the current node
        descNodes = dfsearch( G, numnodes(G), ...
            {'discovernode', 'edgetonew' } );
        
        % Extract edges defining lineage of the current node
        descEdges = descNodes.Edge;
        descEdges( any(isnan(descEdges),2), : ) = [];
        
        descSources = descEdges(:,1);
        descSinks = descEdges(:,2);
        
        % Determine which edge source are division events
        isDiv = ismember(descSources, divStruct.parentID);
        
        % Update segment labels for descendants
        G.Nodes(descSinks,:).Segment = ...
            repmat(newLabel, numel(descSinks), 1);
        
        % Update generation labels for descendants
        for j = 1:size(descEdges, 1)
            
            if isDiv(j)
                G.Nodes(descSinks(j),:).Generation = ...
                    G.Nodes(descSources(j),:).Generation + 1;
            else
                G.Nodes(descSinks(j),:).Generation = ...
                    G.Nodes(descSources(j),:).Generation;
            end
            
        end
        
        % Re-draw display images ------------------------------------------
        T = [prev_T, cur_T, -1];
        
        draw_tracks( G, T, axesList, ImFileName, markerSize);
        
        set(axesList(2), 'XLim', XLim);
        set(axesList(2), 'YLim', YLim);
        
        scaleFigChildren(fig, markerSize);
        
    end

%--------------------------------------------------------------------------
% UPDATE POINT MATCHING
%--------------------------------------------------------------------------
% This function resets the graph structure by running point matching
% STARTING FROM THE CURRENT TIME! All time points up to and including the
% current time will remain unchanged, but all subsequent times will have
% their tracks replaced by simple point matching.

    function updatePointMatch( ~, ~, ~ )
        
        % Confirm that user wishes to run point matching ------------------
        qString = [ 'Are you sure you want to run point matching? ' ...
            'All time points up to and including the current time ' ...
            'will remain unchanged, but all subsequent times will ' ...
            'have their tracks replaced!' ];
        
        answer = questdlg(qString, 'Point Matching Dialog', ...
            'Yes', 'No', 'No' );
        
        if strcmp( answer, 'No' )
            
            userPrompt.String = 'Point matching cancelled';
            
            return;
            
        end
        
        %------------------------------------------------------------------
        % Remove all edges whose target exists at any later time than the
        % current time point
        %------------------------------------------------------------------
        
        cur_T = str2double(timeEdit.String); % The current (child) time
        cur_tidx = find(timePoints == cur_T); % The ID of that time point
        
        if cur_tidx == numel(timePoints)
            
            userPrompt.String = [ 'The current time point is the ' ...
                'final time point. Point matching cannot be updated' ];
            
            return;
            
        end
        
        eID_remove = G.Edges.EndNodes(:,2);
        eID_remove = G.Nodes(eID_remove,:).T;
        eID_remove = find( eID_remove > cur_T );
        
        G = rmedge(G, eID_remove);
        
        %------------------------------------------------------------------
        % Run point matching on all subsequent time points
        %------------------------------------------------------------------
        
        userPrompt.String = 'Point matching 0% complete';
        drawnow;
        
        for tidx = (cur_tidx+1):numel(timePoints)
            
            t = timePoints(tidx);
            
            % Find all of the cell centers at the current time
            curID = find(G.Nodes.T == t);
            curU = G.Nodes(curID,:).UPix;
            
            % Find all of the cell centers at the previous time
            prevID = find(G.Nodes.T == timePoints(tidx-1));
            prevU = G.Nodes(prevID,:).UPix;
            
            % Find the ID of the nearest cell at the previous time to
            % each cell at the current time
            parentID = zeros(size(curID));
            
            for i = 1:numel(curID)
                
                nearID = prevU - repmat( curU(i,:), numel(prevID), 1 );
                nearID = sqrt( sum( nearID.^2, 2 ) );
                [~, nearID] = min(nearID);
                
                parentID(i) = nearID;
                
            end
            
            % Map cell IDs to row IDs in the digraph node table
            parentID = prevID(parentID);
            
            % Update the edge list in the digraph
            G = addedge(G, parentID, curID);
            
            if (mod(tidx, 10) == 0)
                userPrompt.String = ...
                    sprintf( 'Point matching %d%% complete', ...
                    round( 100 * (tidx-(cur_tidx+1)) / ...
                    (numel(timePoints)-(cur_tidx+1)) ) );
                drawnow;
            end
            

        end
        
        userPrompt.String = 'Point matching complete';
        drawnow;
        
        % Update division information
        divStruct = assemble_div_struct(G);
        
        %------------------------------------------------------------------
        % Propagate labels on all subsequent time points
        %------------------------------------------------------------------
        
        % Extract edge endpoints for all times
        sources = G.Edges.EndNodes(:,1);
        sinks = G.Edges.EndNodes(:,2);
        
        % Extract the time point containing each sink
        sinkTime = G.Nodes(sinks,:).T;
        
        userPrompt.String = 'Label propagation 0% complete';
        drawnow;
        
        for tidx = (cur_tidx+1):numel(timePoints)
            
            t = timePoints(tidx);
            
            % Extract the extant nodes at the current time
            nodeIDx = find(G.Nodes.T == t);
            
            % Extract the edges at the current time
            curEdges = (sinkTime == t);
            curSources = sources(curEdges);
            curSinks = sinks(curEdges);
            
            % Set progenitor node generation to 1
            progNodes = ~ismember(nodeIDx, sinks);
            progNodes = nodeIDx(progNodes);
            G.Nodes(progNodes,:).Generation = ones(size(progNodes));
            
            % Determine which edge sources are division events
            isDiv = ismember(curSources, divStruct.parentID);
            
            % Update generation labels
            G.Nodes(curSinks(~isDiv),:).Generation = ...
                G.Nodes(curSources(~isDiv),:).Generation;
            
            G.Nodes(curSinks(isDiv),:).Generation = ...
                G.Nodes(curSources(isDiv),:).Generation + 1;
            
            % Update segment labels
            G.Nodes(curSinks,:).Segment = G.Nodes(curSources,:).Segment;
            
            if (mod(tidx, 10) == 0)
                userPrompt.String = ...
                    sprintf( 'Label propagation %d%% complete', ...
                    round( 100 * (tidx-(cur_tidx+1)) / ...
                    (numel(timePoints)-(cur_tidx+1)) ) );
                drawnow;
            end

        end
        
        userPrompt.String = 'Label propagation complete';
        
        %------------------------------------------------------------------
        % Re-draw display images
        %------------------------------------------------------------------
        T = [-1, cur_T, -1];
        if (cur_tidx > 1), T(1) = timePoints(cur_tidx-1); end
        if (cur_tidx < numel(timePoints))
            T(3) = timePoints(cur_tidx+1);
        end
        
        if (modeMenu.Value == 1)
            
            draw_tracks( G, T, axesList, ImFileName, markerSize);
            
            set(axesList(2), 'XLim', XLim);
            set(axesList(2), 'YLim', YLim);
            
            scaleFigChildren(fig, markerSize);
            
        else
            
            draw_seg( G, T, segLabels, axesList, ImFileName );
            
            set(axesList(2), 'XLim', XLim);
            set(axesList(2), 'YLim', YLim);
            
        end

    end

%--------------------------------------------------------------------------
% DEMONS POINT MATCHING
%--------------------------------------------------------------------------
% Update the lineages for the current parent/child time point combination
% using a Demons-augmented point matching method.
%
% NOTE: THIS METHOD IS ONLY VALID IF ALL IMAGES HAVE IDENTICAL DIMENSIONS
% AND IF THE (x,y)-COORDINATES OF THE PIXEL CENTERS ARE ALREADY IN PIXEL
% SPACE

    function demonsTracking(~, ~, ~)
        
        userPrompt.String = ...
            'This functionality is not yet supported!';
        
        return;
        
        % Confirm that user wishes to run point matching ------------------
        qString = [ 'Are you sure you want to run point matching? ' ...
            'All tracks for the current time point will be replaced.' ];
        
        answer = questdlg(qString, 'Demons Matching Dialog', ...
            'Yes', 'No', 'No' );
        
        if strcmp( answer, 'No' )
            
            userPrompt.String = 'Demons matching cancelled';
            
            return;
            
        end
        
        userPrompt.String = 'Performing Demons matching... ';
        drawnow();
        
        % Extract data about the current/previous time points -------------
        
        cur_T = str2double(timeEdit.String); % The current (child) time
        tidx = find(timePoints == cur_T); % The ID of that time point
        prev_T = timePoints(tidx-1); % The previous (parent) time
        
        % The node IDs of nodes at the current (child) time
        curIDx = find( G.Nodes.T == cur_T );
        
        % The (x,y)-coordinates of all child nodes
        curXY = [ G.Nodes(curIDx, :).X, G.Nodes(curIDx, :).Y ];
        
        % The node IDS of the nodes at the previous (parent) time
        prevIDx = find( G.Nodes.T == prev_T );
        
        % The (x,y)-coordinates of all parent nodes
        prevXY = [ G.Nodes(prevIDx, :).X, G.Nodes(prevIDx, :).Y ];
        
        XLim = get(gca, 'XLim'); % The x-axis limits
        YLim = get(gca, 'YLim'); % The y-axis limits
        
        % Perform Demons-registration on background images ----------------
        
        % Load background images from file
        moving = imadjust( mat2gray( imread( sprintf( ImFileName, ...
            curT ) ) ) );
        target = imadjust( mat2gray( imread( sprintf( ImFileName, ...
            prevT ) ) ) );
        
        imRows = size(moving, 1);
        imCols = size(moving, 2);
        
        % Calculate displacement
        [D, ~] = imregdemons( moving, target, 'DisplayWaitbar', false );
        
        % The (x,y)-displacement fields.  The displacement fields map a
        % pixel in the target image grid to a point in the moving image
        Dx = D(:,:,1);
        Dy = D(:,:,2);
        
        % Update the locations of the nodes at the previous (parent) time
        % using the displacement fields -----------------------------------
        
        % Map cell centers to the nearest pixel center
        regPrevXY = round(prevXY);
        
        regPrevXY( regPrevXY < 1 ) = 1;
        regPrevXY( regPrevXY(:,1) > imRows, 1 ) = imRows;
        regPrevXY( regPrevXY(:,2) > imCols, 2 ) = imCols;
        
        regPrevInd = sub2ind( size(Dx), regPrevXY(:,2), regPrevXY(:,1) );
        
        regPrevXY = [ regPrevXY(:,1) + Dx( regPrevInd ), ...
            regPrevXY(:,2) + Dy( regPrevInd ) ];
        
        % Remove all incoming edges for the nodes at the child time -------
        
        eID_remove = G.Edges.EndNodes(:,2);
        eID_remove = G.Nodes(eID_remove,:).T;
        eID_remove = find( eID_remove == cur_T );
        
        G = rmedge(G, eID_remove);
        
        % Perform point matching on registered locations ------------------
        
        % Find the ID of the nearest cell at the previous time to each cell
        % at the current time
        parentIDx = zeros(size(curIDx));
        
        for i = 1:numel(curIDx)
            
            nearIDx = regPrevXY - repmat( curXY(i,:), numel(prevIDx), 1 );
            nearIDx = sqrt( sum( nearIDx.^2, 2 ) );
            [~, nearIDx] = min(nearIDx);
            
            parentIDx(i) = nearIDx;
            
        end
        
        % Map cell IDs to row IDs in the digraph node table
        parentIDx = prevIDx(parentIDx);
        
        % Update the edge list in the digraph
        G = addedge(G, parentIDx, curIDx);
        
        userPrompt.String = 'Performing Demons matching... Done!';
        
        % Re-draw display images ------------------------------------------
        draw_tracks( G, cur_T, prev_T, prevAx, curAx, fusedAx, ...
            ImFileName, markerSize);
        
        scaleFigChildren(fig, markerSize, XLim, YLim);
        
    end

%==========================================================================
% SEGMENT LABEL MODE FUNTIONS
%==========================================================================

%--------------------------------------------------------------------------
% REMOVE LABEL
%--------------------------------------------------------------------------
% This function removes the currently selected parasegment label in the
% label box and resets all nodes with that label

    function removeLabel( ~, ~, ~ )
        
        if (modeMenu.Value == 1)
            
            userPrompt.String = ['The "Remove Label" function cannot ' ...
                'be used in tracking mode!'];
            
            return;
            
        end
        
        % Determine which label to remove ---------------------------------
        rmID = labelBox.Value;
        rmLabel = segLabels{rmID};
        
        % Confirm that the user wishes to remove the label ----------------
        qString = [ 'Are you sure you want to remove "' rmLabel ...
            '" from the label list? All cells with this label ' ...
            'will be reset!' ];
        
        answer = questdlg( qString, 'Remove Label Dialog', ...
            'Yes', 'No', 'No' );
        
        if strcmp( answer, 'No' )
            userPrompt.String = 'Remove label cancelled';
            return;
        end
        
        % Update the label list -------------------------------------------
        segLabels(rmID) = [];
        
        % Update the label box --------------------------------------------
        
        % Determine which nodes should contribute to the coloring scheme
        hasLabel = ismember(G.Nodes.Segment, segLabels);
        
        % The number of unique generations
        numGen = max(G.Nodes(hasLabel,:).Generation);
        
        % Assign a unique color to each label/generation
        segColors = 255 .* distinguishable_colors( ...
            numel(segLabels)+numGen, [0 0 0; 1 1 1] );
        
        genColors = segColors(1:numGen, :);
        segColors = segColors((numGen+1):end, :);
        
        % Set the generation string
        genTypes = { 'PSPR', 'Post-Wave 1', ...
            'Post-Wave 2', 'Diff Cleavage %d'};
        genString = cell(size(genColors,1),1);
        for i = 1:size(genColors,1)
            
            if (i < 4)
                genString{i} = sprintf( ...
                    ['<HTML><BODY bgcolor = "rgb(%f,%f,%f">' ...
                    genTypes{i} '</BODY></HTML>'], ...
                    genColors(i,1), genColors(i,2), genColors(i,3) );
            else
                genString{i} = sprintf( ...
                    ['<HTML><BODY bgcolor = "rgb(%f,%f,%f">' ...
                    genTypes{4} '</BODY></HTML>'], ...
                    genColors(i,1), genColors(i,2), genColors(i,3), i-3 );
            end
            
        end
        
        genBox.String = genString;
        
        % Set the label string
        labelString = cell(size(segLabels));
        for i = 1:numel(segLabels)
            labelString{i} = sprintf( ...
                ['<HTML><BODY bgcolor="rgb(%f,%f,%f)">' ...
                segLabels{i} '</BODY></HTML>'], ...
                segColors(i,1), segColors(i,2), segColors(i,3) );
        end
        
        labelBox.String = labelString;
        
        if isempty(segLabels)
            labelBox.Value = [];
        else
            labelBox.Value = 1;
        end
        
        % Update the graph structure --------------------------------------
        rmIDx = strcmp( rmLabel, G.Nodes.Segment );
        
        if any(rmIDx)
            G.Nodes(rmIDx,:).Segment = repmat({'none'}, sum(rmIDx), 1);
        end
        
        % Re-draw display images ------------------------------------------
        
        % The ID of the time point to display
        tidx = timeSld.Value;
        
        % The times of interest
        T = [ -1, timePoints(tidx), -1];
        if (tidx > 1), T(1) = timePoints(tidx-1); end
        if (tidx < numel(timePoints)), T(3) = timePoints(tidx+1); end
        
        % Get current axis limits
        XLim = get(axesList(2), 'XLim');
        YLim = get(axesList(2), 'YLim');
        
        % Render images
        draw_seg( G, T, segLabels, axesList, ImFileName );
        
        % Set axes limits
        set(axesList(2), 'XLim', XLim);
        set(axesList(2), 'YLim', YLim);
        
        userPrompt.String = ['"' rmLabel '" removed from the label list'];
        
    end

%--------------------------------------------------------------------------
% ADD LABEL
%--------------------------------------------------------------------------
% Add a label to the label list

    function addLabel( ~, ~, ~ )
        
        if (modeMenu.Value == 1)
            
            userPrompt.String = ['The "Add Label" function cannot ' ...
                'be used in tracking mode!'];
            
            return;
            
        end
        
        % Prompt for user input -------------------------------------------
        
        prompt = {'Enter new label:'};
        dlgtitle = 'New Label Input';
        dims = [1 35];
        
        newLabel = inputdlg( prompt, dlgtitle, dims );
        if isempty( newLabel ), return; end
        newLabel = newLabel{1};
        
        % Update label list -----------------------------------------------
        
        % Check if the label is already on the list
        if any( strcmp( newLabel, segLabels ) )
            
            userPrompt.String = [ '"' newLabel '" is already a ' ...
                'member of the label list' ];
            
            return;
            
        end
        
        segLabels{end+1} = newLabel;
        
        % Update the label box --------------------------------------------
        
        % Determine which nodes should contribute to the coloring scheme
        hasLabel = ismember(G.Nodes.Segment, segLabels);
        
        % The number of unique generations
        if ~any(hasLabel)
            numGen = 1;
        else
            numGen = max(G.Nodes(hasLabel,:).Generation);
        end
        
        % Assign a unique color to each label/generation
        segColors = 255 .* distinguishable_colors( ...
            numel(segLabels)+numGen, [0 0 0; 1 1 1] );
        
        genColors = segColors(1:numGen, :);
        segColors = segColors((numGen+1):end, :);
        
        % Set the generation string
        genTypes = { 'PSPR', 'Post-Wave 1', ...
            'Post-Wave 2', 'Diff Cleavage %d'};
        genString = cell(size(genColors,1),1);
        for i = 1:size(genColors,1)
            
            if (i < 4)
                genString{i} = sprintf( ...
                    ['<HTML><BODY bgcolor = "rgb(%f,%f,%f">' ...
                    genTypes{i} '</BODY></HTML>'], ...
                    genColors(i,1), genColors(i,2), genColors(i,3) );
            else
                genString{i} = sprintf( ...
                    ['<HTML><BODY bgcolor = "rgb(%f,%f,%f">' ...
                    genTypes{4} '</BODY></HTML>'], ...
                    genColors(i,1), genColors(i,2), genColors(i,3), i-3 );
            end
            
        end
        
        genBox.String = genString;
        
        % Set the label string
        labelString = cell(size(segLabels));
        for i = 1:numel(segLabels)
            labelString{i} = sprintf( ...
                ['<HTML><BODY bgcolor="rgb(%f,%f,%f)">' ...
                segLabels{i} '</BODY></HTML>'], ...
                segColors(i,1), segColors(i,2), segColors(i,3) );
        end
        
        labelBox.String = labelString;
        
        userPrompt.String = ['"' newLabel '" added to label list'];

    end

%--------------------------------------------------------------------------
% SET LABEL
%--------------------------------------------------------------------------
% Sets a label for a node at the current time and propagates that label to
% all of that nodes descendants

    function setLabel( ~, ~, ~ )
        
        if (modeMenu.Value == 1)
            
            userPrompt.String = ['The "Set Label" function cannot ' ...
                'be used in tracking mode!'];
            
            return;
            
        end
        
        % Determine if the user wants to overwrite existing labels --------
        qString = [ 'Do you wish to overwrite the existing labels ' ...
            'of the desecendants of the selected nodes?' ];
        
        answer = questdlg( qString, 'Set Label Dialog', ...
            'Yes', 'No', 'No' );
        
        if strcmp( answer, 'Yes' )
            overwriteLabel = true;
        else
            overwriteLabel = false;
        end
        
        % Extract data about the current time points ----------------------
        
        % The ID of the time point being displayed
        tidx = timeSld.Value;
        
        % The times of interest
        T = [ -1, timePoints(tidx), -1];
        if (tidx > 1), T(1) = timePoints(tidx-1); end
        if (tidx < numel(timePoints)), T(3) = timePoints(tidx+1); end
        
        % The node IDs of the nodes at the current time
        curIDx = find( G.Nodes.T == T(2) );
        
        % The locations of the current nodes in pixel space
        curU = G.Nodes(curIDx, :).UPix;
        
        % Get the patch graphics object displayed at the current time
        curPatches = get(axesList(2), 'Children');
        curPatches = curPatches(strcmp(get(curPatches, 'Type'), 'patch'));
        
        % The patch that displays the segment label
        curSegPatch = curPatches(2);
        
        XLim = get(axesList(2), 'XLim'); % The x-axis limits
        YLim = get(axesList(2), 'YLim'); % The y-axis limits
        
        % Extract the active label
        labelID = labelBox.Value;
        activeLabel = segLabels{labelID};
        
        % Determine which nodes should contribute to the coloring scheme
        hasLabel = ismember(G.Nodes.Segment, segLabels);
        
        % The number of unique generations
        if any(hasLabel)
            numGen = max(G.Nodes(hasLabel,:).Generation);
        else
            numGen = 0;
        end
        
        % Assign a unique color to the active label
        labelColor = distinguishable_colors( ...
            numel(segLabels)+numGen, [0 0 0; 1 1 1] );
        labelColor = labelColor(numGen+labelID, :);
        
        % Modify user interaction strings
        userPrompt.String = [ 'Please select cells to give label "' ...
            activeLabel '". Push button to indicate when ' ...
            'selection is complete' ];
        
        if overwriteLabel
            userPrompt.String = [ userPrompt.String ...
                '. Overwriting labels' ];
        else
            userPrompt.String = [ userPrompt.String ...
                '. Not overwriting labels' ];
        end
        
        userActionBtn.String = 'Finish Selection';
        
        % Re-set the user action button in case the user forgot
        userActionBtn.Value = 0;
        
        % Extract cells to endow with the active label --------------------
        
        % The node IDs that will receive the label
        addIDx = [];
        curAddIDx = [];
        
        count = 0; % An indexing variable
        while true
            
            % The (x,y)-coordinates of the current selection
            [xC, yC] = ginput(1);
            drawnow();
            
            % Check the the specified point lies within the axis limits
            if ( (xC <= XLim(1)) || (XLim(2) <= xC) || ...
                    (yC <= YLim(1)) || (YLim(2) <= yC) )
 
                if userActionBtn.Value
                    break;
                else
                    userPrompt.String = [ 'Please choose a point that '...
                        'lies within the axes limits!' ];
                    continue;
                end

            end
            
            % Point match to find the node ID of the current selection
            curAddID = curU - repmat([xC, yC], size(curU,1), 1);
            curAddID = sqrt( sum( curAddID.^2, 2 ) );
            [~, curAddID] = min(curAddID);
            addID = curIDx(curAddID);
            
            if ismember(addID, addIDx)
                userPrompt.String = [ 'Selected node has already ' ...
                    'been picked. Please choose another' ];
                continue;
            end
            
            % Increment the index variable
            count = count + 1;
            
            % Add the current selection to the list of newly labelled nodes
            addIDx(count) = addID;
            curAddIDx(count) = curAddID;
            
            % Visualize the current selection
            curSegPatch.FaceVertexCData(curAddID, :) = labelColor;
            curSegPatch.FaceVertexAlphaData(curAddID, :) = 0.2;
            
        end
        
        % Re-set the user action button
        userActionBtn.Value = 0;
        
        % Re-set visualization for clarity
        curSegPatch.FaceVertexCData(curAddIDx, :) = ...
            0.8 .* ones( numel(addIDx), 3 );
        curSegPatch.FaceVertexAlphaData(curAddIDx, :) = 0;
        
        % Update interaction strings
        userActionBtn.String = [];
       
        if ~isempty(addIDx)
            userPrompt.String = 'Selection finished. Updating labels... ';
        else
            userPrompt.String = 'No nodes selected';
            return;
        end
        
        % Update the graph structure --------------------------------------
        
        for i = 1:numel(addIDx)
            
           % Extract all descendants of the current node
           descNodes = dfsearch( G, addIDx(i), ...
               {'discovernode', 'edgetonew' } );
           descNodes = descNodes.Node;
           descNodes = descNodes( ~isnan(descNodes) );
           
           G.Nodes(descNodes(1),:).Segment = {activeLabel};
           
           for j = 2:numel(descNodes)
               
               if overwriteLabel
                   G.Nodes(descNodes(j),:).Segment = {activeLabel};
               elseif strcmp(G.Nodes(descNodes(j),:).Segment, 'none')
                   G.Nodes(descNodes(j),:).Segment = {activeLabel};
               end
               
           end
           
           userPrompt.String = ...
                sprintf( [ 'Selection finished. Updating labels... ' ...
                '%d%% complete' ], round(100 * i / numel(addIDx)) );
            
            drawnow;
           
        end
        
        % Update label boxes ----------------------------------------------
        
        % Determine which nodes should contribute to the coloring scheme
        hasLabel = ismember(G.Nodes.Segment, segLabels);
        
        % The number of unique generations
        if any(hasLabel)
            numGen = max(G.Nodes(hasLabel,:).Generation);
        else
            numGen = 0;
        end
        
        % Assign a unique color to each label/generation
        segColors = 255 .* distinguishable_colors( ...
            numel(segLabels)+numGen, [0 0 0; 1 1 1] );
        
        genColors = segColors(1:numGen, :);
        segColors = segColors((numGen+1):end, :);
        
        % Set the generation string
        genTypes = { 'PSPR', 'Post-Wave 1', ...
            'Post-Wave 2', 'Diff Cleavage %d'};
        genString = cell(size(genColors,1),1);
        for i = 1:size(genColors,1)
            
            if (i < 4)
                genString{i} = sprintf( ...
                    ['<HTML><BODY bgcolor = "rgb(%f,%f,%f">' ...
                    genTypes{i} '</BODY></HTML>'], ...
                    genColors(i,1), genColors(i,2), genColors(i,3) );
            else
                genString{i} = sprintf( ...
                    ['<HTML><BODY bgcolor = "rgb(%f,%f,%f">' ...
                    genTypes{4} '</BODY></HTML>'], ...
                    genColors(i,1), genColors(i,2), genColors(i,3), i-3 );
            end
            
        end
        
        genBox.String = genString;
        
        % Set the label string
        labelString = cell(size(segLabels));
        for i = 1:numel(segLabels)
            labelString{i} = sprintf( ...
                ['<HTML><BODY bgcolor="rgb(%f,%f,%f)">' ...
                segLabels{i} '</BODY></HTML>'], ...
                segColors(i,1), segColors(i,2), segColors(i,3) );
        end
        
        labelBox.String = labelString;
        
        % Re-draw display images ------------------------------------------
        
        % Render images
        draw_seg( G, T, segLabels, axesList, ImFileName );
        
        % Set axes limits
        set(axesList(2), 'XLim', XLim);
        set(axesList(2), 'YLim', YLim);
        
        userPrompt.String = 'Selection finished. Updating labels... Done';
        
    end

%--------------------------------------------------------------------------
% UNSET LABEL
%--------------------------------------------------------------------------
% Remove the label from a user selected noded and all of its descendants

    function unsetLabel( ~, ~, ~ )
        
        if (modeMenu.Value == 1)
            
            userPrompt.String = ['The "Unset Label" function cannot ' ...
                'be used in tracking mode!'];
            
            return;
            
        end
        
        % Determine if the user wants to overwrite existing labels --------
        qString = [ 'Do you wish to remove the existing labels ' ...
            'of the desecendants of the selected node?' ];
        
        answer = questdlg( qString, 'Unset Label Dialog', ...
            'Yes', 'No', 'No' );
        
        if strcmp( answer, 'Yes' )
            overwriteLabel = true;
        else
            overwriteLabel = false;
        end
        
        % Extract data about the current time points ----------------------
        
        % The ID of the time point being displayed
        tidx = timeSld.Value;
        
        % The times of interest
        T = [ -1, timePoints(tidx), -1];
        if (tidx > 1), T(1) = timePoints(tidx-1); end
        if (tidx < numel(timePoints)), T(3) = timePoints(tidx+1); end
        
        % The node IDs of the nodes at the current time
        curIDx = find( G.Nodes.T == T(2) );
        
        % The locations of the current nodes in pixel space
        curU = G.Nodes(curIDx, :).UPix;
        
        XLim = get(axesList(2), 'XLim'); % The x-axis limits
        YLim = get(axesList(2), 'YLim'); % The y-axis limits
        
        % Select the node for removal -------------------------------------
        
        % Prompt for user input
        userPrompt.String = [ 'Please select the cell whose ' ...
            'label will be removed' ];
        
        if overwriteLabel
            userPrompt.String = [ userPrompt.String ...
                '. Removing label of descendants' ];
        else
            userPrompt.String = [ userPrompt.String ...
                '. Not removing labels of descendants' ];
        end
        
        [xC, yC] = ginput(1);
        
        % Check the the specified point lies within the axis limits
        if ( (xC <= XLim(1)) || (XLim(2) <= xC) || ...
                (yC <= YLim(1)) || (YLim(2) <= yC) )
            
            userPrompt.String = [ 'Please choose a point that '...
                'lies within the axes limits!' ];
            
            return;
            
        end
        
        % Point match to find the correct node ID
        rmID = curU - repmat( [xC, yC], size(curU,1), 1);
        rmID = sqrt( sum( rmID.^2, 2 ) );
        [~, rmID] = min(rmID);
        rmID = curIDx(rmID);
        
        % Update the graph structure --------------------------------------
        
        userPrompt.String = 'Selection finished. Updating labels...';
        
        % Extract all descendants of the current node
        descNodes = dfsearch( G, rmID, ...
            {'discovernode', 'edgetonew' } );
        descNodes = descNodes.Node;
        descNodes = descNodes( ~isnan(descNodes) );
        
        G.Nodes(descNodes(1),:).Segment = {'none'};
        
        for j = 2:numel(descNodes)
            if overwriteLabel
                G.Nodes(descNodes(j),:).Segment = {'none'};
            end
        end
        
        % Update label boxes ----------------------------------------------
        
        % Determine which nodes should contribute to the coloring scheme
        hasLabel = ismember(G.Nodes.Segment, segLabels);
        
        % The number of unique generations
        numGen = max(G.Nodes(hasLabel,:).Generation);
        
        % Assign a unique color to each label/generation
        segColors = 255 .* distinguishable_colors( ...
            numel(segLabels)+numGen, [0 0 0; 1 1 1] );
        
        genColors = segColors(1:numGen, :);
        segColors = segColors((numGen+1):end, :);
        
        % Set the generation string
        genTypes = { 'PSPR', 'Post-Wave 1', ...
            'Post-Wave 2', 'Diff Cleavage %d'};
        genString = cell(size(genColors,1),1);
        for i = 1:size(genColors,1)
            
            if (i < 4)
                genString{i} = sprintf( ...
                    ['<HTML><BODY bgcolor = "rgb(%f,%f,%f">' ...
                    genTypes{i} '</BODY></HTML>'], ...
                    genColors(i,1), genColors(i,2), genColors(i,3) );
            else
                genString{i} = sprintf( ...
                    ['<HTML><BODY bgcolor = "rgb(%f,%f,%f">' ...
                    genTypes{4} '</BODY></HTML>'], ...
                    genColors(i,1), genColors(i,2), genColors(i,3), i-3 );
            end
            
        end
        
        genBox.String = genString;
        
        % Set the label string
        labelString = cell(size(segLabels));
        for i = 1:numel(segLabels)
            labelString{i} = sprintf( ...
                ['<HTML><BODY bgcolor="rgb(%f,%f,%f)">' ...
                segLabels{i} '</BODY></HTML>'], ...
                segColors(i,1), segColors(i,2), segColors(i,3) );
        end
        
        labelBox.String = labelString;
        
        % Re-draw display images ------------------------------------------
        
        % Render images
        draw_seg( G, T, segLabels, axesList, ImFileName );
        
        % Set axes limits
        set(axesList(2), 'XLim', XLim);
        set(axesList(2), 'YLim', YLim);
        
        userPrompt.String = 'Selection finished. Updating labels... Done';
        
    end

%--------------------------------------------------------------------------
% PROPAGATE GENERATION LABELS
%--------------------------------------------------------------------------
% This function resets the generation labels by propagating generations
% down lineages as derived from the tracking structure.  If you had a
% perfect tracking, a single run of this function would produce a perfect
% result.

    function propGen( ~, ~, ~ )
        
        if (modeMenu.Value == 1)
            
            userPrompt.String = ['The "Propagate Generations" ' ...
                'function cannot be used in tracking mode!'];
            
            return;
            
        end
        
        % Confirm that the user wishes to propagate generation labels -----
        
        qString = [ 'Are you sure that you want to propagate '...
            'generations labels?  All time points up to and ' ...
            'including the current time will reman unchanged, ' ...
            'but all subsequent times will have their generation ' ...
            'labels replaced!' ];
        
        answer = questdlg( qString, 'Propagate Generation Dialog', ...
            'Yes', 'No', 'No' );
        
        if strcmp( answer, 'No' )
            
            userPrompt.String = 'Generation propagation cancelled';
            
            return;
            
        end
        
        % Extract date about the current time points ----------------------
        
        % The ID of the time point being displayed
        tidx = timeSld.Value;
        
        % The times of interest
        T = [ -1, timePoints(tidx), -1];
        if (tidx > 1), T(1) = timePoints(tidx-1); end
        if (tidx < numel(timePoints)), T(3) = timePoints(tidx+1); end
        
        % Extract current time point axis limits
        XLim = get(axesList(2), 'XLim'); % The x-axis limits
        YLim = get(axesList(2), 'YLim'); % The y-axis limits
        
        % Propagate generation labels on all subsequent time points -------
        
        % Extract edge endpoints for all times
        sources = G.Edges.EndNodes(:,1);
        sinks = G.Edges.EndNodes(:,2);
        
        % Extract the time point containing each sink
        sinkTime = G.Nodes(sinks,:).T;
        
        userPrompt.String = 'Generation propagation 0% complete';
        drawnow;
        
        for ttidx = (tidx+1):numel(timePoints)

            t = timePoints(ttidx);
            
            % Extract the extant nodes at the current time
            nodeIDx = find(G.Nodes.T == t);
            
            % Extract the edges at the current time
            curEdges = (sinkTime == t);
            curSources = sources(curEdges);
            curSinks = sinks(curEdges);
            
            % Set progenitor node generation to 1
            progNodes = ~ismember(nodeIDx, sinks);
            progNodes = nodeIDx(progNodes);
            G.Nodes(progNodes,:).Generation = ones(size(progNodes));
            
            % Determine which edge sources are division events
            isDiv = ismember(curSources, divStruct.parentID);
            
            % Update generation labels
            G.Nodes(curSinks(~isDiv),:).Generation = ...
                G.Nodes(curSources(~isDiv),:).Generation;
            
            G.Nodes(curSinks(isDiv),:).Generation = ...
                G.Nodes(curSources(isDiv),:).Generation + 1;
            
            if (mod(ttidx, 10) == 0)
                userPrompt.String = ...
                    sprintf( 'Generation propagation %d%% complete', ...
                    round( 100 * (ttidx-(tidx+1)) / ...
                    (numel(timePoints)-(tidx+1)) ) );
                drawnow;
            end
            
        end
        
        userPrompt.String = 'Generation propagation complete';
        
        % Update label boxes ----------------------------------------------
        
        % Determine which nodes should contribute to the coloring scheme
        hasLabel = ismember(G.Nodes.Segment, segLabels);
        
        % The number of unique generations
        numGen = max(G.Nodes(hasLabel,:).Generation);
        
        % Assign a unique color to each label/generation
        segColors = 255 .* distinguishable_colors( ...
            numel(segLabels)+numGen, [0 0 0; 1 1 1] );
        
        genColors = segColors(1:numGen, :);
        segColors = segColors((numGen+1):end, :);
        
        % Set the generation string
        genTypes = { 'PSPR', 'Post-Wave 1', ...
            'Post-Wave 2', 'Diff Cleavage %d'};
        genString = cell(size(genColors,1),1);
        for i = 1:size(genColors,1)
            
            if (i < 4)
                genString{i} = sprintf( ...
                    ['<HTML><BODY bgcolor = "rgb(%f,%f,%f">' ...
                    genTypes{i} '</BODY></HTML>'], ...
                    genColors(i,1), genColors(i,2), genColors(i,3) );
            else
                genString{i} = sprintf( ...
                    ['<HTML><BODY bgcolor = "rgb(%f,%f,%f">' ...
                    genTypes{4} '</BODY></HTML>'], ...
                    genColors(i,1), genColors(i,2), genColors(i,3), i-3 );
            end
            
        end
        
        genBox.String = genString;
        
        % Set the label string
        labelString = cell(size(segLabels));
        for i = 1:numel(segLabels)
            labelString{i} = sprintf( ...
                ['<HTML><BODY bgcolor="rgb(%f,%f,%f)">' ...
                segLabels{i} '</BODY></HTML>'], ...
                segColors(i,1), segColors(i,2), segColors(i,3) );
        end
        
        labelBox.String = labelString;
        
        % Re-draw display images ------------------------------------------
        
        % Render images
        draw_seg( G, T, segLabels, axesList, ImFileName );
        
        % Set axes limits
        set(axesList(2), 'XLim', XLim);
        set(axesList(2), 'YLim', YLim);
        
    end

%--------------------------------------------------------------------------
% SET GENERATION LABELS
%--------------------------------------------------------------------------
% Set the generation label for a node at the current tiem and propagate
% that label to all of the nodes descendants accounting for division events

    function setGen( ~, ~, ~ )
        
        if (modeMenu.Value == 1)
            
            userPrompt.String = ['The "Set Generations" ' ...
                'function cannot be used in tracking mode!'];
            
            return;
            
        end
        
        % Determine if the user wants to overwrite existing labels --------
        qString = [ 'Do you wish to overwrite the existing generation' ...
            ' labels of the descendants of the selected nodes?' ];
        
        answer = questdlg( qString, 'Set Generation Dialog', ...
            'Yes', 'No', 'No' );
        
        if strcmp( answer, 'Yes' )
            overwriteLabel = true;
        else
            overwriteLabel = false;
        end
        
        % Determine the new generation label to assign to selected points -
        pString = ['Please provide the generation label that ' ...
            'will be assigned to selected nodes.'];
        
        newGen = inputdlg( pString, 'Set Generation Dialog', ...
            [1 35], {'-1'});
        
        newGen = str2double(newGen);
        
        if ((newGen < 1) || (mod(newGen,1) ~= 0))
            
            userPrompt.String = ['New generation label must be a ' ...
                'positive integer >= 1'];
            
            return;
            
        end
        
        % Extract data about the current time points ----------------------
        
        % The ID of the time point being displayed
        tidx = timeSld.Value;
        
        % The times of interest
        T = [ -1, timePoints(tidx), -1];
        if (tidx > 1), T(1) = timePoints(tidx-1); end
        if (tidx < numel(timePoints)), T(3) = timePoints(tidx+1); end
        
        % The node IDs of the nodes at the current time
        curIDx = find( G.Nodes.T == T(2) );
        
        % The locations of the current nodes in pixel space
        curU = G.Nodes(curIDx, :).UPix;
        
        % Get the patch graphics object displayed at the current time
        curPatches = get(axesList(2), 'Children');
        curPatches = curPatches(strcmp(get(curPatches, 'Type'), 'patch'));
        
        % The patch that displays the generation label
        % curGenPatch = curPatches(1); 
        
        % The patch that displays the segment label
        curSegPatch = curPatches(2);
        
        % Extract current time axis limits
        XLim = get(axesList(2), 'XLim'); % The x-axis limits
        YLim = get(axesList(2), 'YLim'); % The y-axis limits
        
        % Modify user interaction strings
        userPrompt.String = [ 'Please select cells to give ' ...
            'generation "' num2str(newGen) '". ' ...
            'Push button to indicate when selection is complete' ];
        
        if overwriteLabel
            userPrompt.String = [ userPrompt.String ...
                '. Overwriting labels' ];
        else
            userPrompt.String = [ userPrompt.String ...
                '. Not overwriting labels' ];
        end
        
        userActionBtn.String = 'Finish Selection';
        
        % Re-set the user action button in case the user forgot
        userActionBtn.Value = 0;
        
        drawnow;
        
        % Extract cells to endow with the new generation label ------------
        
        % The node IDs that will receive the label
        addIDx = [];
        curAddIDx = [];
        curAddColors = [];
        
        count = 0; % An indexing variable
        while true
            
            % The (x,y)-coordinates of the current selection
            [xC, yC] = ginput(1);
            drawnow();
            
            % Check the the specified point lies within the axis limits
            if ( (xC <= XLim(1)) || (XLim(2) <= xC) || ...
                    (yC <= YLim(1)) || (YLim(2) <= yC) )
 
                if userActionBtn.Value
                    break;
                else
                    userPrompt.String = [ 'Please choose a point that '...
                        'lies within the axes limits!' ];
                    continue;
                end

            end
            
            % Point match to find the node ID of the current selection
            curAddID = curU - repmat([xC, yC], size(curU,1), 1);
            curAddID = sqrt( sum( curAddID.^2, 2 ) );
            [~, curAddID] = min(curAddID);
            addID = curIDx(curAddID);
            
            if ismember(addID, addIDx)
                userPrompt.String = [ 'Selected node has already ' ...
                    'been picked. Please choose another' ];
                continue;
            end
            
            % Increment the index variable
            count = count + 1;
            
            % Add the current selection to the list of newly labelled nodes
            addIDx(count) = addID;
            curAddIDx(count) = curAddID;
            curAddColors(count,:) = ...
                curSegPatch.FaceVertexCData(curAddID, :);
            
            % Visualize the current selection
            curSegPatch.FaceVertexCData(curAddID, :) = [0 0 0];
            
        end
        
        % Re-set the user action button
        userActionBtn.Value = 0;
        
        % Update interaction strings
        userActionBtn.String = [];
       
        if ~isempty(addIDx)
            userPrompt.String = 'Selection finished. Updating labels... ';
            drawnow;
        else
            userPrompt.String = 'No nodes selected';
            return;
        end
        
        % Update the graph structure --------------------------------------
        
        for i = 1:numel(addIDx)

            G.Nodes(addIDx(i),:).Generation = newGen;
            
            if overwriteLabel
                
                % Extract all descendants of the current node
                descNodes = dfsearch(G, addIDx(i), ...
                    {'discovernode', 'edgetonew'} );
                
                % Extract edges defining lineage of the current node
                descEdges = descNodes.Edge;
                descEdges( any(isnan(descEdges),2), : ) = [];
                
                descSources = descEdges(:,1);
                descSinks = descEdges(:,2);
                
                % Determine which edge source are division events
                isDiv = ismember(descSources, divStruct.parentID);
                
                % Update generation labels for descendants
                for j = 1:size(descEdges, 1)
                    
                    if isDiv(j)
                        G.Nodes(descSinks(j),:).Generation = ...
                            G.Nodes(descSources(j),:).Generation + 1;
                    else
                        G.Nodes(descSinks(j),:).Generation = ...
                            G.Nodes(descSources(j),:).Generation;
                    end
                    
                end

            end
            
            userPrompt.String = ...
                sprintf( [ 'Selection finished. Updating labels... ' ...
                '%d%% complete' ], round(100 * i / numel(addIDx)) );
            
            drawnow;

        end
        
        % Update label boxes ----------------------------------------------
        
        % Determine which nodes should contribute to the coloring scheme
        hasLabel = ismember(G.Nodes.Segment, segLabels);
        
        % The number of unique generations
        numGen = max(G.Nodes(hasLabel,:).Generation);
        
        % Assign a unique color to each label/generation
        segColors = 255 .* distinguishable_colors( ...
            numel(segLabels)+numGen, [0 0 0; 1 1 1] );
        
        genColors = segColors(1:numGen, :);
        segColors = segColors((numGen+1):end, :);
        
        % Set the generation string
        genTypes = { 'PSPR', 'Post-Wave 1', ...
            'Post-Wave 2', 'Diff Cleavage %d'};
        genString = cell(size(genColors,1),1);
        for i = 1:size(genColors,1)
            
            if (i < 4)
                genString{i} = sprintf( ...
                    ['<HTML><BODY bgcolor = "rgb(%f,%f,%f">' ...
                    genTypes{i} '</BODY></HTML>'], ...
                    genColors(i,1), genColors(i,2), genColors(i,3) );
            else
                genString{i} = sprintf( ...
                    ['<HTML><BODY bgcolor = "rgb(%f,%f,%f">' ...
                    genTypes{4} '</BODY></HTML>'], ...
                    genColors(i,1), genColors(i,2), genColors(i,3), i-3 );
            end
            
        end
        
        genBox.String = genString;
        
        % Set the label string
        labelString = cell(size(segLabels));
        for i = 1:numel(segLabels)
            labelString{i} = sprintf( ...
                ['<HTML><BODY bgcolor="rgb(%f,%f,%f)">' ...
                segLabels{i} '</BODY></HTML>'], ...
                segColors(i,1), segColors(i,2), segColors(i,3) );
        end
        
        labelBox.String = labelString;
        
        % Re-draw display images ------------------------------------------
        
        % Render images
        draw_seg( G, T, segLabels, axesList, ImFileName );
        
        % Set axes limits
        set(axesList(2), 'XLim', XLim);
        set(axesList(2), 'YLim', YLim);
        
        userPrompt.String = 'Selection finished. Updating labels... Done';
        
    end


%==========================================================================
% GENERAL FUNCTIONS
%==========================================================================

%--------------------------------------------------------------------------
% KEY PRESS FUNCTION
%--------------------------------------------------------------------------
% A catch-all function for key presses while the GUI window is active.
%   - Shift the display time using the arrow keys

    function kp_Fcn( ~, event )
        
        % The key that was pressed
       keyPress = event.Key;
       
       % The list of keys that will trigger an action
       activeKeys = {'leftarrow', 'uparrow', 'downarrow', ...
           'rightarrow', 'w', 'a', 's', 'd'};
       
       % Do nothing if an active key was not pressed
       if ~ismember( keyPress, activeKeys ), return; end
       
       switch keyPress
           
           % Switch to segment label mode ---------------------------------
           case {'downarrow', 's'}
               
               % Get the current time point of interest
               tidx = timeSld.Value;
               
               % Update time points of interest
               T = [ -1, timePoints(tidx), -1];
               if (tidx > 1), T(1) = timePoints(tidx-1); end
               if (tidx < numel(timePoints))
                   T(3) = timePoints(tidx+1);
               end
               
               if (modeMenu.Value == 2)
                   
                   return;
                   
               else
                   
                   modeMenu.Value = 2;
                   
                   % Update image labels
                   axLabel1.String = sprintf('T = %d', T(1));
                   axLabel2.String = sprintf('T = %d', T(2));
                   axLabel3.String = sprintf('T = %d', T(3));
                   
                   % Get current axis limits
                   XLim = get(axesList(2), 'XLim');
                   YLim = get(axesList(2), 'YLim');
                   
                   % Render images
                   draw_seg( G, T, segLabels, axesList, ImFileName );
                   
                   % Set axes limits
                   set(axesList(2), 'XLim', XLim);
                   set(axesList(2), 'YLim', YLim);
                   
               end
               
               return;
               
           % Switch to tracking mode --------------------------------------
           case {'w', 'uparrow'}
               
               % Get the current time point of interest
               tidx = timeSld.Value;
               
               % Update time points of interest
               T = [ -1, timePoints(tidx), -1];
               if (tidx > 1), T(1) = timePoints(tidx-1); end
               if (tidx < numel(timePoints))
                   T(3) = timePoints(tidx+1);
               end
               
               if ((modeMenu.Value == 1) || (tidx == 1))
                   
                   return;
                   
               else
                   
                   modeMenu.Value = 1;
                   
                   % Update image labels
                   axLabel1.String = sprintf('T = %d', T(1));
                   axLabel2.String = 'FUSED';
                   axLabel3.String = sprintf('T = %d', T(2));
                   
                   % Get current axis limits
                   XLim = get(axesList(2), 'XLim');
                   YLim = get(axesList(2), 'YLim');
                   
                   % Render images
                   draw_tracks(G, T, axesList, ImFileName, markerSize)
                   
                   % Set axes limits
                   set(axesList(2), 'XLim', XLim);
                   set(axesList(2), 'YLim', YLim);
                   
                   % Re-scale figure children
                   scaleFigChildren(fig, markerSize);
                   
               end
               
               return;
               
           % Increment the display time -----------------------------------
           case {'d', 'rightarrow'}
               
               % The ID of the current display time
               tidx = timeSld.Value;
               
               % Do nothing if the increment would go out of bounds
               if tidx == numel(timePoints), return; end
               
               % Increment the index
               tidx = tidx+1;
               
               % Update slider value
               timeSld.Value = tidx;
               
               % The times of interest
               T = [ -1, timePoints(tidx), -1];
               if (tidx > 1), T(1) = timePoints(tidx-1); end
               if (tidx < numel(timePoints))
                   T(3) = timePoints(tidx+1);
               end
               
               % Update the edit control display
               timeEdit.String = num2str(T(2));
               
               if (modeMenu.Value == 1)
                   
                   % Update image labels
                   axLabel1.String = sprintf('T = %d', T(1));
                   axLabel2.String = ' ';
                   axLabel3.String = sprintf('T = %d', T(2));
                   
                   % Get current axis limits
                   XLim = get(axesList(2), 'XLim');
                   YLim = get(axesList(2), 'YLim');
                   
                   % Render images
                   draw_tracks(G, T, axesList, ImFileName, markerSize)
                   
                   % Set axes limits
                   set(axesList(2), 'XLim', XLim);
                   set(axesList(2), 'YLim', YLim);
                   
                   % Re-scale figure children
                   scaleFigChildren(fig, markerSize);
                   
               else
                   
                   % Update image labels
                   axLabel1.String = sprintf('T = %d', T(1));
                   axLabel2.String = sprintf('T = %d', T(2));
                   axLabel3.String = sprintf('T = %d', T(3));
                   
                   % Get current axis limits
                   XLim = get(axesList(2), 'XLim');
                   YLim = get(axesList(2), 'YLim');
                   
                   % Render images
                   draw_seg( G, T, segLabels, axesList, ImFileName );
                   
                   % Set axes limits
                   set(axesList(2), 'XLim', XLim);
                   set(axesList(2), 'YLim', YLim);
                   
               end
               
               return;
               
           % Decrement the display time -----------------------------------
           case {'a', 'leftarrow'}
               
               % The ID of the current display time
               tidx = timeSld.Value;
               
               % Do nothing if the deccrement would go out of bounds
               if tidx == 1, return; end
               if ((modeMenu.Value == 1) && (tidx == 2)), return; end
               
               % Decrement the index
               tidx = tidx-1;
               
               % Update slider value
               timeSld.Value = tidx;
               
               % The times of interest
               T = [ -1, timePoints(tidx), -1];
               if (tidx > 1), T(1) = timePoints(tidx-1); end
               if (tidx < numel(timePoints)) 
                   T(3) = timePoints(tidx+1);
               end
               
               % Update the edit control display
               timeEdit.String = num2str(T(2));
               
               if (modeMenu.Value == 1)
                   
                   % Update image labels
                   axLabel1.String = sprintf('T = %d', T(1));
                   axLabel2.String = ' ';
                   axLabel3.String = sprintf('T = %d', T(2));
                   
                   % Get current axis limits
                   XLim = get(axesList(2), 'XLim');
                   YLim = get(axesList(2), 'YLim');
                   
                   % Render images
                   draw_tracks(G, T, axesList, ImFileName, markerSize)
                   
                   % Set axes limits
                   set(axesList(2), 'XLim', XLim);
                   set(axesList(2), 'YLim', YLim);
                   
                   % Re-scale figure children
                   scaleFigChildren(fig, markerSize);
                   
               else
                   
                   % Update image labels
                   axLabel1.String = sprintf('T = %d', T(1));
                   axLabel2.String = sprintf('T = %d', T(2));
                   axLabel3.String = sprintf('T = %d', T(3));
                   
                   % Get current axis limits
                   XLim = get(axesList(2), 'XLim');
                   YLim = get(axesList(2), 'YLim');
                   
                   % Render images
                   draw_seg( G, T, segLabels, axesList, ImFileName );
                   
                   % Set axes limits
                   set(axesList(2), 'XLim', XLim);
                   set(axesList(2), 'YLim', YLim);
                   
               end
               
               return;
               
           % Do nothing if no active keys were pressed --------------------
           otherwise
               return;
               
       end
        
    end

%--------------------------------------------------------------------------
% RESET ZOOM
%--------------------------------------------------------------------------
% This function resets the the axis limits of each figure to their full
% extent and restores the function of the 'double-click' zoom operation -
% essentially it does what the 'double-click' zoom operation is supposed to
% do automatically.  It should be noted that this callback function is only
% necessary because MATLAB is absolute trash.

    function resetZoom( ~, ~, ~ )
        
        % The the current time axis
        ax = axesList(2);
        
        % Get the children of the current axis
        axKids = get(ax, 'Children');
        
        % Extract the size of the background image
        for i = 1:length(axKids)
            if strcmp( get(axKids(i), 'Type'), 'image' )
                
                XLim = get(axKids(i), 'XData');
                YLim = get(axKids(i), 'YData');
                
                break;
            end
        end
        
        % Reset the zoom to the size of the image
        set(axesList(2), 'XLim', XLim);
        set(axesList(2), 'YLim', YLim);
        
        % Re-scale figure children if the GUI is in tracking mode
        if (modeMenu.Value == 1)
            scaleFigChildren(fig, markerSize);
        end
        
        zoom reset
        
    end

%--------------------------------------------------------------------------
% ZOOM CALLBACK
%--------------------------------------------------------------------------
% When any zoom operation is performed, this function makes sure
% the marker size of any scatter plots are changed to reflect the scaling

    function zoomCallBack( ~, ~ )
        
        scaleFigChildren(fig, markerSize);

    end

%--------------------------------------------------------------------------
% UPDATE TIME SLIDER
%--------------------------------------------------------------------------
% This function updates the current times displayed in the GUI based on
% activation of the slider control

    function updateTimeSld( ~, ~, ~ )
        
        % The ID of the time point to display
        % ( Match to nearest valid time point index )
        tidx = timeSld.Value;
        
        % Point match the input to the nearest valid time point index
        allTIDx = (1:numel(timePoints)).';
        tidx = allTIDx - repmat(tidx, numel(allTIDx), 1);
        tidx = sqrt( sum( tidx.^2, 2 ) );
        [~, tidx] = min(tidx);
        tidx = allTIDx(tidx);
        
        if ((modeMenu.Value == 1) && (tidx == 1))
            
            userPrompt.String = ['The first time point cannot be ' ...
                'treated as a child time point in tracking mode. ' ...
                'Using second time point instead' ];
            
            tidx = 2;
            
        end
        
        % Update slider value
        timeSld.Value = tidx;
        
        % The times of interest
        T = [ -1, timePoints(tidx), -1];
        if (tidx > 1), T(1) = timePoints(tidx-1); end
        if (tidx < numel(timePoints)), T(3) = timePoints(tidx+1); end
        
        % Update the edit control display
        timeEdit.String = num2str(T(2));
        
        if (modeMenu.Value == 1)
            
            % Update image labels
            axLabel1.String = sprintf('T = %d', T(1));
            axLabel2.String = ' ';
            axLabel3.String = sprintf('T = %d', T(2));
            
            % Get current axis limits
            XLim = get(axesList(2), 'XLim');
            YLim = get(axesList(2), 'YLim');
            
            % Render images
            draw_tracks(G, T, axesList, ImFileName, markerSize)
            
            % Set axes limits
            set(axesList(2), 'XLim', XLim);
            set(axesList(2), 'YLim', YLim);
            
            % Re-scale figure children
            scaleFigChildren(fig, markerSize);
            
        else
            
            % Update image labels
            axLabel1.String = sprintf('T = %d', T(1));
            axLabel2.String = sprintf('T = %d', T(2));
            axLabel3.String = sprintf('T = %d', T(3));
            
            % Get current axis limits
            XLim = get(axesList(2), 'XLim');
            YLim = get(axesList(2), 'YLim');
            
            % Render images
            draw_seg( G, T, segLabels, axesList, ImFileName );
            
            % Set axes limits
            set(axesList(2), 'XLim', XLim);
            set(axesList(2), 'YLim', YLim);
            
        end
        
    end

%--------------------------------------------------------------------------
% UPDATE TIME EDIT
%--------------------------------------------------------------------------
% This function updates the current times displayed in the GUI based on the
% activation of the edit control

    function updateTimeEdit( ~, ~, ~ )
        
        % The times of interest
        T = [-1, str2double(timeEdit.String), -1];
        
        % Point match to find the nearest valid time point index
        tidx = timePoints -repmat( T(2), numel(timePoints), 1);
        tidx = sqrt( sum( tidx.^2, 2 ) );
        [~, tidx] = min(tidx);
        
        if ((modeMenu.Value == 1) && (tidx == 1))
            
            userPrompt.String = ['The first time point cannot be ' ...
                'treated as a child time point in tracking mode. ' ...
                'Using second time point instead' ];
            
            tidx = 2;
            
        end
        
        % Update time points of interest
        T = [ -1, timePoints(tidx), -1];
        if (tidx > 1), T(1) = timePoints(tidx-1); end
        if (tidx < numel(timePoints)), T(3) = timePoints(tidx+1); end
        
        % Update edit value
        timeEdit.String = num2str(T(2));
        
        % Update the slider control
        timeSld.Value = tidx;
        
        if (modeMenu.Value == 1)
            
            % Update image labels
            axLabel1.String = sprintf('T = %d', T(1));
            axLabel2.String = 'FUSED';
            axLabel3.String = sprintf('T = %d', T(2));
            
            % Get current axis limits
            XLim = get(axesList(2), 'XLim');
            YLim = get(axesList(2), 'YLim');
            
            % Render images
            draw_tracks(G, T, axesList, ImFileName, markerSize)
            
            % Set axes limits
            set(axesList(2), 'XLim', XLim);
            set(axesList(2), 'YLim', YLim);
            
            % Re-scale figure children
            scaleFigChildren(fig, markerSize);
            
        else
            
            % Update image labels
            axLabel1.String = sprintf('T = %d', T(1));
            axLabel2.String = sprintf('T = %d', T(2));
            axLabel3.String = sprintf('T = %d', T(3));
            
            % Get current axis limits
            XLim = get(axesList(2), 'XLim');
            YLim = get(axesList(2), 'YLim');
            
            % Render images
            draw_seg( G, T, segLabels, axesList, ImFileName );
            
            % Set axes limits
            set(axesList(2), 'XLim', XLim);
            set(axesList(2), 'YLim', YLim);
            
        end
        
    end

%--------------------------------------------------------------------------
% GUI MODE SELECTION
%--------------------------------------------------------------------------
% This function switches the GUI between tracking mode and segment label
% mode

    function setGUIMode( ~, ~, ~ )
        
        % Get the current time point of interest
        tidx = timeSld.Value;
        
        % Update time points of interest
        T = [ -1, timePoints(tidx), -1];
        if (tidx > 1), T(1) = timePoints(tidx-1); end
        if (tidx < numel(timePoints)), T(3) = timePoints(tidx+1); end
        
        if ((modeMenu.Value == 1) && (tidx > 1))
            
            % Update image labels
            axLabel1.String = sprintf('T = %d', T(1));
            axLabel2.String = 'FUSED';
            axLabel3.String = sprintf('T = %d', T(2));
            
            % Get current axis limits
            XLim = get(axesList(2), 'XLim');
            YLim = get(axesList(2), 'YLim');
            
            % Render images
            draw_tracks(G, T, axesList, ImFileName, markerSize)
            
            % Set axes limits
            set(axesList(2), 'XLim', XLim);
            set(axesList(2), 'YLim', YLim);
            
            % Re-scale figure children
            scaleFigChildren(fig, markerSize);
            
        else
            
            % Update image labels
            axLabel1.String = sprintf('T = %d', T(1));
            axLabel2.String = sprintf('T = %d', T(2));
            axLabel3.String = sprintf('T = %d', T(3));
            
            % Get current axis limits
            XLim = get(axesList(2), 'XLim');
            YLim = get(axesList(2), 'YLim');
            
            % Render images
            draw_seg( G, T, segLabels, axesList, ImFileName );
            
            % Set axes limits
            set(axesList(2), 'XLim', XLim);
            set(axesList(2), 'YLim', YLim);
            
        end
        
    end

%--------------------------------------------------------------------------
% EXIT GENERATION GUIDE
%--------------------------------------------------------------------------
% This function checks for consistency in the output graph structure and
% closes the picker GUI

    function exitGuide( ~, ~, ~ )
        
        % Display output message
        userPrompt.String = 'Goodbye!';
        
        %*****************************************
        % DO CONSISTENCY CHECKS - WHAT ARE THEY???
        pause(1); % I want them to see the message
        %*****************************************
        
        close all force
        
    end

%--------------------------------------------------------------------------
% SAVE PROGRESS
%--------------------------------------------------------------------------
% Save progress without exiting the GUI

    function saveProgress( ~, ~, ~ )
       
        userPrompt.String = 'Saving progress... ';
        save('Parhyale_Master_GUI_Progress', 'G');
        userPrompt.String = 'Saving progress... Done';
        
    end


end

%**************************************************************************
% DRAW TRACKS
%**************************************************************************
function draw_tracks(G, T, axesList, ImFileName, markerSize)
%DRAW_TRACKS Draws the figures depicting the tracks and lineages
%
%   INPUT PARAMETERS:
%
%       - G:        The updated MATLAB digraph object representating the
%                   manually curated tracks and lineages
%
%       - T:        The three time points currently under consideration:
%                   [ previousT currentT nextT]
%
%       - axesList:     A graphics array containing the axes used to
%                       display the images: { previous current next }
%
%       - ImFileName:   The name of the image file over which to
%                       display the parasegments at each given
%                       time point.  It is assumed that each time
%                       point is saved separately with identical
%                       names up to an integer which delineated the
%                       time point.
%
%       - markerSize:   The length scale of each marker face in units
%                       based on the x-axis of each displayed image

%--------------------------------------------------------------------------
% Validate Inputs
%--------------------------------------------------------------------------

% Extract time points
timePoints = unique(G.Nodes.T);

if ~ismember(T(2), timePoints(2:end))
    warning('Invalid child time point supplied for rendering');
end

%--------------------------------------------------------------------------
% Load/Extract Data
%--------------------------------------------------------------------------

% Load background images from file
curMIP = imadjust( mat2gray( imread( sprintf( ImFileName, T(2) ) ) ) );
prevMIP = imadjust( mat2gray( imread( sprintf( ImFileName, T(1) ) ) ) );

% The extant nodes/locations at the child time point
curID = find(G.Nodes.T == T(2));
curU = G.Nodes(curID,:).UPix;

% The extant nodes/locations at the parent time point
prevID = find(G.Nodes.T == T(1));
prevU = G.Nodes(prevID,:).UPix;

% The end nodes of each edge in the tracking graph
sources = G.Edges.EndNodes(:,1);
sinks = G.Edges.EndNodes(:,2);

% Logical vector indicating if each current node has a parent
[isSink, edgeLoc] = ismember(curID, sinks);

% Extract the parent node of each current node (There should only be one!)
parentID = zeros(size(curID));
parentID(isSink) = sources(edgeLoc(isSink));

% Map parent node IDs onto indices into the reduced node list
parentID = changem( parentID, (0:numel(prevID)).', [0; prevID] );

% Find the number of nodes at the current time with no parent node
numProg = sum(parentID == 0);

% Add a new ID for each progenitor node at the current time
parentID( parentID == 0 ) = numel(prevID) + (1:numProg).';

%--------------------------------------------------------------------------
% Generate the Figures
%--------------------------------------------------------------------------

% Establish a separate color for each distinct lineage (i.e. the number of
% progenitor cells)
linColors = distinguishable_colors( numel(prevID)+numProg, ...
    [0 0 0; 1 1 1; 0 1 0; 1 0 1] );

% Generate the previous (parent) time plot --------------------------------
subplot(axesList(1));

imshow(prevMIP);
hold on
scatter( prevU(:,1), prevU(:,2), [], ...
    linColors(1:numel(prevID), :), 's' );
hold off
set(gca, 'YDir', 'normal');

% Generate the fused image plot -------------------------------------------
subplot(axesList(2));

imshowpair(curMIP, prevMIP);
hold on
scatter( prevU(:,1), prevU(:,2), [], ...
    linColors(1:numel(prevID), :), 's' );
scatter( curU(:,1), curU(:,2), [], ...
    linColors(parentID, :), 'x' );
hold off
set(gca, 'YDir', 'normal');

% Generate the current (child) time plot ----------------------------------
subplot(axesList(3));

imshow(curMIP);
hold on
scatter( curU(:,1), curU(:,2), [], ...
    linColors(parentID, :), 'x' );
hold off
set(gca, 'YDir', 'normal');

% Re-scale the markers in all figure children -----------------------------
fig = get(axesList(2), 'Parent');
scaleFigChildren(fig, markerSize);

end

%**************************************************************************
% SCALE FIGURE CHILDREN
%**************************************************************************
function scaleFigChildren(fig, markerSize)
%SCALEFIGCHILDREN This function re-scales scatter plot markers to fit the
%limits of the display axes
%
%   INPUT PARAMETERS:
%
%       - fig:          The handle to the GUI figure window
%
%       - markerSize:   The linear size of the scatter plot marker
%

% Get the current axes limits
XLim = get(gca, 'XLim');
% YLim = get(gca, 'YLim');

% Get the children of the current figure window
figKids = get(fig, 'Children');

for k = 1:length(figKids)
    if strcmp( get(figKids(k), 'Type'), 'axes' )
        
        % Set the current axis
        ax = figKids(k);
        
        % Set the current axes units to points
        currentunits = get(ax, 'Units');
        set(ax, 'Units', 'Points');
        
        % Get the position of the current axes in point
        axPos = get(ax, 'Position');
        
        % Reset axis units
        set(ax, 'Units', currentunits);
        
        % Calculate marker width
        markerWidth = markerSize/diff(XLim)*axPos(3);
        
        % Get the children of the current axes
        axKids = get(ax, 'Children');
        
        % Re-scale all scatter plot markers
        for i = 1:length(axKids)
            if strcmp( get(axKids(i), 'Type'), 'scatter' )
                set( axKids(i), 'SizeData', markerWidth.^2 );
                set( axKids(i), 'LineWidth', markerWidth/8 );
            end
        end
        
    end
end

end

%**************************************************************************
% DRAW SEGMENTS
%**************************************************************************
function draw_seg( G, T, segLabels, axesList, ImFileName )
%DRAW_SEG Draws the figures depicting the parasegments
%
%   INPUT PARAMETERS:
%
%       - G:        The updated MATLAB digraph object representating the
%                   manually curated tracks and lineages
%
%       - T:        The three time points currently under consideration:
%                   [ previousT currentT nextT]
%
%       - segLabels:    A cell array containing the unique parasegment
%                       labels
%
%       - axesList:     A graphics array containing the axes used to
%                       display the images: { previous current next }
%
%       - ImFileName:   The name of the image file over which to
%                       display the parasegments at each given
%                       time point.  It is assumed that each time
%                       point is saved separately with identical
%                       names up to an integer which delineated the
%                       time point.

%--------------------------------------------------------------------------
% Validate Inputs
%--------------------------------------------------------------------------

% Extract time points
timePoints = unique(G.Nodes.T);

if ~ismember(T(2), timePoints)
    warning('Invalid time points supplied for rendering');
    return;
end

%--------------------------------------------------------------------------
% Extract Parasegment Label/Cell Cycle Generation Lists
%--------------------------------------------------------------------------

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

%--------------------------------------------------------------------------
% Iterate Over Axes to Create Images
%--------------------------------------------------------------------------

for i = 1:3
    
    t = T(i); % The current time point being rendered
    subplot(axesList(i)); % The current axis for rendering
    
    %----------------------------------------------------------------------
    % Handle Out-of-Bounds Neighbor Time Points
    %----------------------------------------------------------------------
    if ~ismember( t, timePoints )
        
        MIP = imread( sprintf( ImFileName, timePoints(1) ) );
        MIP = zeros(size(MIP));
        
        imshow(MIP);
        
        continue;
        
    end
    
    %----------------------------------------------------------------------
    % Load/Extract Data
    %----------------------------------------------------------------------
    
    % Load background image from file
    MIP = imadjust( mat2gray( imread( sprintf( ImFileName, t ) ) ) );
    
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
    
    maxFaceSize = max(cellfun(@(x) numel(x), c));
    voronoiFace = nan(size(c,1), maxFaceSize);
    for k = 1:size(c,1)
        voronoiFace(k, 1:numel(c{k})) = c{k};
    end
    
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
            fSegAlpha( locNode ) = 1;
            dispNodes( locNode ) = true;
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
    
    % Basic background image
    imshow(MIP);
    
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
    
end

end

%**************************************************************************
% ASSEMBLE DIVISION EVENT INFORMATION
%**************************************************************************
function divStruct = assemble_div_struct(G, rmBadDiv)
% ASSEMBLE_DIV_STRUCT Detects and conglomerates relative information about
% cell divisions divined from the input tracking structure
%
%   INPUT PARAMETERS:
%
%       - G:        The updated MATLAB digraph object representating the
%                   manually curated tracks and lineages
%
%       - rmBadDiv  If true, division events not meeting particular quality
%                   checks are removed from consideration
%
%   OUTPUT PARAMETERS:
%
%       - divStruct:    The structure containing the division data

if (nargin < 2), rmBadDiv = true; end

sources = G.Edges.EndNodes(:,1);
sinks = G.Edges.EndNodes(:,2);

dupl = struct2cell(find_duplicate_rows(sources)).';

% The node IDs of all parent nodes in division events
divParentIDx = dupl(:,1);

% The node IDs of the daughter cells in the corresponding division events
divChildIDx = cellfun( @(x) sinks(x), dupl(:,2), 'UniformOutput', false);

% Remove all events with more than two daughter cells (false positives)
badEvent =  cellfun( @numel, divChildIDx ) ~= 2;

if rmBadDiv
    
    divParentIDx(badEvent) = [];
    divChildIDx(badEvent) = [];
    
end

% Extract the time of the division event
divTime = cellfun( @(x) G.Nodes(x,:).T, divParentIDx );

% Convert cells to arrays
divParentIDx = cell2mat(divParentIDx);
divChildIDx = [ divChildIDx{:} ].';

% Extract the division axis angles (relative to image x-axis)
U12 = G.Nodes(divChildIDx(:,2), :).UPix - ...
    G.Nodes(divChildIDx(:,1), :).UPix;
U12( U12(:,1) < 0, : ) =  -U12( U12(:,1) < 0, : );
U12 = U12 ./ sqrt( sum( U12.^2, 2 ) );

divAngle = atan2(U12(:,2), U12(:,1));

% Determine division event labels
divParentLabel = G.Nodes(divParentIDx,:).Segment;
divChildLabel = [ G.Nodes(divChildIDx(:,1),:).Segment, ...
    G.Nodes(divChildIDx(:,2),:).Segment ];

% Remove events where both children do have the same label as their parent
if rmBadDiv
    
    badEvent = strcmp( divChildLabel(:,1), divChildLabel(:,2) );
    badEvent = badEvent & strcmp( divChildLabel(:,1), divParentLabel );
    badEvent = ~badEvent;
    
    divParentIDx(badEvent) = [];
    divChildIDx(badEvent,:) = [];
    divTime(badEvent) = [];
    divAngle(badEvent) = [];
    divParentLabel(badEvent) = [];
    
end

% Construct a struct that contains all the division information
divStruct = struct( 'parentID', divParentIDx, 'childID', divChildIDx, ...
    't', divTime, 'angle', divAngle, 'label', {divParentLabel} );

end

%**************************************************************************
% FIND DUPLICATE ROWS
%**************************************************************************
function [ dupl, C ] = find_duplicate_rows( A )
%FIND_DUPLICATE_ROWS Detect duplicate rows in a numeric matrix
%
%   INPUT PARAMETERS:
%
%       - A         MxN numeric matrix
%
%   OUTPUT PARAMETERS:
%
%       - dupl      A structure containing the value of the duplicate rows
%                   and the locations of the duplicate rows in the original
%                   input matrix
%
%       - C         M'xN' numeric matrix containing the sorted unique rows
%                   contained in A

[C, ia, ~] = unique( A, 'rows' );

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
    dupl_idx{i} = find( abs(sum( diffA, 2 )) < eps );
    
end

dupl = struct( 'val', dupl_val, 'idx', dupl_idx );

end
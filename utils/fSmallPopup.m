function iIndRefImg = fSmallPopup(cShow,sTitle,iMax,iSelected)
% show popup figure for an easy selection (RegGUI, EvalGUI)
%
% input:
% cShow             cell array with strings for shown listbox
% sTitle            figure title
% iMax              amount of to be selectable items
%
% output:
% iIndRefImg        selected items from shown listbox
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner
% -------------------------------------------------------------------------

iFIGUREWIDTH = 300;
iFIGUREHEIGHT = 400;
iBUTTONHEIGHT = 24;

iIndRefImg = [];
iPos = get(0, 'ScreenSize');
if(nargin < 4), iSelected = 1; end;
if(nargin < 3), iMax = 1; end;

% -------------------------------------------------------------------------
% Create figure and GUI elements
hF = figure( ...
    'Position'              , [(iPos(3) - iFIGUREWIDTH)/2, (iPos(4) - iFIGUREHEIGHT)/2, iFIGUREWIDTH, iFIGUREHEIGHT], ...
    'Units'                 , 'pixels', ...
    'DockControls'          , 'off', ...
    'WindowStyle'           , 'modal', ...
    'Name'                  , sTitle, ...
    'NumberTitle'           , 'off', ...
    'Resize'                , 'off');

hList = uicontrol(hF, ...
    'Style'                 , 'listbox', ...
    'Units'                 , 'pixels', ...
    'Position'              , [1 iBUTTONHEIGHT + 1 iFIGUREWIDTH iFIGUREHEIGHT - iBUTTONHEIGHT], ...
    'ForegroundColor'       , 'w' , ...
    'BackgroundColor'       , 'k', ...
    'HitTest'               , 'on', ...
    'String'                , cShow, ...
    'Min'                   , 0, ...
    'Max'                   , iMax, ...
    'Value'                 , iSelected, ...
    'Callback'              , @fMouseActionFcn);

hButOK = uicontrol(hF, ...
    'Style'                 , 'pushbutton', ...
    'Units'                 , 'pixels', ...
    'Position'              , [1 1 iFIGUREWIDTH/2 iBUTTONHEIGHT], ...
    'ForegroundColor'       , 'w', ...
    'BackgroundColor'       , 'k', ...
    'Callback'              , @fMouseActionFcn, ...
    'HitTest'               , 'on', ...
    'String'                , 'OK');

hButCancel = uicontrol(hF, ...
    'Style'                 , 'pushbutton', ...
    'Units'                 , 'pixels', ...
    'Position'              , [iFIGUREWIDTH/2 + 1 1 iFIGUREWIDTH/2 iBUTTONHEIGHT], ...
    'ForegroundColor'       , 'w', ...
    'BackgroundColor'       , 'k', ...
    'Callback'              , @fMouseActionFcn, ...
    'HitTest'               , 'on', ...
    'String'                , 'Cancel');
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Set default action and enable gui interaction
sAction = 'Cancel';
uiwait(hF);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% uiresume was triggered (in fMouseActionFcn) -> return
try
    if strcmp(sAction, 'OK')
        iIndRefImg = get(hList, 'Value');
%         csVarOut = cell(length(iList), 1);
%         for iI = 1:length(iList)
%             csVarOut(iI) = csVars(iList(iI));
%         end
    end
    close(hF);
catch %#ok<CTCH>
    iIndRefImg = [];
end
% -------------------------------------------------------------------------


    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * *
    % * * NESTED FUNCTION fMouseActionFcn (nested in fGetRefImg)
    % * *
    % * * Determine whether axes are linked
	% * *
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    function fMouseActionFcn(hObject, eventdata)
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % React on action depending on its source component
        switch(hObject)
            
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Click in LISTBOX: return if double-clicked
            case hList
                if strcmp(get(hF, 'SelectionType'), 'open')
                    sAction = 'OK';
                    uiresume(hF);
                end
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % OK button
            case hButOK
                sAction = 'OK';
                uiresume(hF);
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % CANCEL button
            case hButCancel
                sAction = 'Cancel';
                uiresume(hF);
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            otherwise

        end
        % End of switch statement
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % * * END NESTED FUNCTION fGridMouseMoveFcn
	% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    end
end
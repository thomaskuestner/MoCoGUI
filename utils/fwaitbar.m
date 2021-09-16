function h = fwaitbar( iVal, h )
% waitbar wrapper to change its appearance
%
% input:
% iVal              waitbar value
% sMessage          waitbar text
%
% output:
% h                 waitbar handle
%
% -------------------------------------------------------------------------
% (c) 2015: Thomas Kuestner
% -------------------------------------------------------------------------

if(nargin < 2)
    h = waitbar(iVal); % just update current waitbar
else 
    if ischar(h) || iscellstr(h) % init
        sMessage = h;
        h = waitbar(iVal,sMessage);
        set(h,'Color','k','Name','Progress');
        haxes = get(findobj(h,'Type','Axes'));
        set(haxes.Title,'Color','w','FontSize',12);
        drawnow expose;
    elseif all(ishghandle(h, 'figure')) % update
        h = waitbar(iVal,h);        
    end
end

end


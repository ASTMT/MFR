% File: vlUtResizeFigure.m
%
% Syntax: nrUtResizeFigure(fighandle,width,height)
%
% Description:
%       Resizes a figure OuterPosition to the width, height specified in
%       pixels.  Useful for resizing figures to ensure that the size is
%       sufficient to not cut off numbers or text in legends.  Adjusts
%       position so that figure grows towards the bottom of the screen.
%
% Input Parameters:
%       fighandle - the handle to the figure to me modified. gcf can be
%       supplied for the currently active figure.
%       width - new figure width in pixels
%       height - new figure height in pixels
%
% Output Parameters:
%       none.
%
% Required Global Data Structures:
%       none.
%
% Required Data Files:
%       none.
%       

% Original version written by Chris Carter (ccarter@tmt.org)
% Tuesday 16th October 2012

%
% Revision History
%
% $Id: vlUtResizeFigure.m,v 1.2 2012/10/19 21:31:45 roberts Exp $
%
% INDENT-OFF*
% $Log: vlUtResizeFigure.m,v $
% Revision 1.2  2012/10/19 21:31:45  roberts
% removed windowstyle modifications
%
% Revision 1.1  2012/10/19 16:09:51  roberts
% Initial revision
%
%
% INDENT-ON*

function vlUtResizeFigure(fighandle,width,height)

% set(0,'DefaultFigureWindowStyle','normal');
% 
% set(fighandle,'WindowStyle','normal');

position = get(gcf,'OuterPosition');

ypos = position(2) - (height - position(4));
if ypos < 100
    ypos = 100;
end

set(fighandle,'OuterPosition',[position(1),ypos,width,height]);

end

% End of file


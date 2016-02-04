% File: <vlSciPlotHist.m>
%
% Syntax: handle = vlSciPlotHist(vec,numbin,plottitle)
%
% Description:
%       Plots a histogram into numbin bins from the input data vector
%
% Input Parameters:
%       vec     - data vector
%       numbin  - number of bins
%       plottitle - title to add to plot
%       
% Output Parameters:
%       handle - figure handle
%
% Required Global Data Structures:
%       none
%
% Required Data Files:
%       none
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% static char rcsid[] = "$Id: vlSciPlotHist.m,v 1.2 2006/01/08 22:55:49 roberts Exp $";
% INDENT-OFF*
%
% INDENT-ON*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%           Herzberg Institute of Astrophysics                  %%%%%
%%%%%%      Astronomy Technology Research Group - Victoria           %%%%%
%
% (c) <2003>				        (c) <2003>
% National Research Council		    Conseil national de recherches
% Ottawa, Canada, K1A 0R6 		    Ottawa, Canada, K1A 0R6
% All rights reserved			    Tous droits reserves
% 					
% NRC disclaims any warranties,	    Le CNRC denie toute garantie
% expressed, implied, or statu-	    enoncee, implicite ou legale,
% tory, of any kind with respect	de quelque nature que se soit,
% to the software, including		concernant le logiciel, y com-
% without limitation any war-		pris sans restriction toute
% ranty of merchantability or		garantie de valeur marchande
% fitness for a particular pur-	    ou de pertinence pour un usage
% pose.  NRC shall not be liable	particulier.  Le CNRC ne
% in any event for any damages,	    pourra en aucun cas etre tenu
% whether direct or indirect,		responsable de tout dommage,
% special or general, consequen-	direct ou indirect, particul-
% tial or incidental, arising		ier ou general, accessoire ou
% from the use of the software.	    fortuit, resultant de l'utili-
% 					                sation du logiciel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handle = vlSciPlotHist(vec,numbin,plottitle)

%ticks = round(numbin/10);
ticks = 1;
maxval = max(vec);
minval = min(vec);

binval = zeros(numbin,1);
cum = zeros(numbin,1);
binamt = zeros(numbin,1);

for ii = 1 : numbin
    binval(ii) = ii*(maxval-minval)/numbin + minval;
    cum(ii) = sum(vec <= binval(ii));
end

for ii = 2 : numbin
    binamt(ii) = cum(ii) - cum(ii-1);
end
binamt(1) = cum(1);

%xaxisval = round(binval*10)/10;
xaxisval = binval;

% Plot histogram
handle = figure;
bar(binamt);

title(sprintf('%s, %d items, [max, min = %.3e,%.3e]',plottitle,length(vec),maxval,minval));
%xlabel('Range of Values');
ylabel('# Values');
set(gca,'XLimMode','manual');
xlim([0 numbin+1]);
set(gca,'XTickMode','manual');
set(gca,'XTick',1:ticks:numbin);
% Set the X axis tick labels by creating a cell array of strings
xvals = [];
for i = 1:ticks:numbin
    if isempty(xvals)
        xvals = sprintf('%+12.3e',xaxisval(i));
    else
        xvals = [xvals; sprintf('%+12.3e',xaxisval(i))];
    end
end

% Format and rotate the X tick labels
% See Matlab Technical Solution Number: 1-15TK6
% How can I rotate my X-axis tick labels and place an X-label on my plot?

ax=axis;    % current axis limits
Yl = ax(3:4); % Y-axis limits
Xt=1:ticks:numbin;
axis(axis); % Set the axis limit modes (e.g. XLimMode) to manual
t = text(Xt,Yl(1)*ones(1,length(Xt)),xvals(:,1:end));
set(t,'HorizontalAlignment','right','VerticalAlignment','top', 'Rotation',90);
% Remove the default labels
set(gca,'XTickLabel','')

%set(gca,'XTickLabel',xvals);
%set(gca,'XTickLabel',xaxisval(1:ticks:numbin));


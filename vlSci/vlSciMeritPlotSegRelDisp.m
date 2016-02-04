% File: <vlSciMeritPlotSegRelDisp.m>
%
% Syntax: function [h,xyuv,max_shear,max_gap,min_gap,gslist] = vlSciMeritPlotSegRelDisp(xyuv,SegSize,PlotScale,plotflag)
%
% Description:
%       Makes a plot of the primary mirror segment relative displacement
%
% Input Parameters:
%       xyuv    -   (n x 4). Position (x,y) and displacement (u,v) of each
%                   segment.
%       SegSize    - (scalar) Segment diameter to determine the nearest
%                   neighbours.
%       PlotScale   - Scale to plot segment displacement vectors.
%       plotflag  - (binary) 1 if figure should be created, 0 if only xyuv
%                   is calculated and returned.
%
% Output Parameters:
%       h   - handle to the figure object
%       uvxy - input uvxy with missing segment rows removed
%       max_shear - maximum relative shear between segment edges
%       max_gap - largest increase in relative segment gap distance
%       min_gap - largest decrease in relative segment gap distance
%       gslist - list of mid point positions between segments, with
%                       gap and shear data (X,Y,gap,shear).

%
% Required Global Data Structures:
%       none.
%
% Required Data Files:
%       none.
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: vlSciMeritPlotSegRelDisp.m,v 1.6 2012/07/24 17:23:22 roberts Exp $
%
% INDENT-OFF*
% $Log: vlSciMeritPlotSegRelDisp.m,v $
% Revision 1.6  2012/07/24 17:23:22  roberts
% added flag to control whether plot is created
%
% Revision 1.5  2006/07/18 20:34:41  roberts
% modified to return data on all segments
%
% Revision 1.4  2006/04/13 22:11:29  roberts
% fixed bug when bset or rset are empty
%
% Revision 1.3  2006/01/17 23:28:02  roberts
% Added figure command to open new figure window
%
% Revision 1.2  2006/01/08 22:51:59  roberts
% updated for R2 release
%
% Revision 1.1  2006/01/08 18:25:43  roberts
% Initial revision
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

function [h,xyuv,max_shear,max_gap,min_gap,gslist] = vlSciMeritPlotSegRelDisp(xyuv,SegSize,PlotScale,plotflag)

% remove zero rows from xyuv
ind = (xyuv == 0);
ind = (sum(ind')' == size(xyuv,2));
xyuv(ind,:) = [];

rset = [];
bset = [];
gap = [];
shear = [];

for seg = 1:length(xyuv) 
    % find positions relative to this segment
    relpos=(xyuv(:,1:2)-repmat(xyuv(seg,1:2),length(xyuv),1));
    % find distances of other segments from this seg
    dist=(relpos(:,1).^2+relpos(:,2).^2).^.5;
    % select only adjacent segments
    adjsegs = xyuv(dist<SegSize&dist>0,:);
    segdist = dist(dist<SegSize&dist>0,:);
    for i = 1:size(adjsegs)
        mag = adjsegs(i,3:4) - xyuv(seg,3:4);
        newdist = adjsegs(i,1:2) - xyuv(seg,1:2) + mag;
        pos = (xyuv(seg,1:2) + adjsegs(i,1:2))/2;
        angle = acos(dot(newdist,mag)/(norm(newdist)*norm(mag)));
        gap = norm(mag)*cos(angle);
        shear = norm(mag)*sin(angle);
        if norm(newdist) > segdist(i) 
            rset = [rset; pos mag gap shear];
        else
            bset = [bset; pos mag gap shear];
        end
    end
end

if plotflag
    figure;
    h = scatter(xyuv(:,1),xyuv(:,2),'g');
    hold on
else
    h = 0;
end

if ~isempty(rset)
    if plotflag
        quiver(rset(:,1),rset(:,2),rset(:,3),rset(:,4),PlotScale,'r');
    end
else
    rset = [0 0 0 0 0 0];
end

if ~isempty(bset)
    if plotflag
        quiver(bset(:,1),bset(:,2),bset(:,3),bset(:,4),PlotScale,'b');
    end
else
    bset = [0 0 0 0 0 0];
end
max_shear = max([max(rset(:,6)) max(bset(:,6))]);
max_gap = max(rset(:,5));
min_gap = min(bset(:,5));

if plotflag
    hold off
end

if isequal(bset,[0 0 0 0 0 0]) && isequal(rset,[0 0 0 0 0 0])
	error(‘both best and rset are empty’);
end

if isequal(bset,[0 0 0 0 0 0])
	gslist = rset;
elseif isequal(rest,[0 0 0 0 0 0])
	gslist = bset;
else
	gslist = [bset ; rset];
end
% remove the mag fields
gslist(:,3:4) = [];
% remove repeated values
gslist = unique(gslist,'rows');

% end of vlSciMeritPlotSegRelDisp.m
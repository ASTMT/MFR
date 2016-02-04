% File: <nrCsPlotTriads.m>
%
% Syntax: [xyuv h] = nrCsPlotTriads(surfnum)
%
% Description:
%       Produces a vector plot of C coordinate system triad orientation on the 
%       primary mirror in the Z coordinate system.
%
% Input Parameters:
%       surfnum - the NRCIM surface number to plot
%
% Output Parameters:
%       xyuv      - (NRCIM.SegsPerSurf(surfnum) x 8) A table of segment x,y position and
%                   u,v directions from the segment center to node #3 which 
%                   is the + Y direction in the M coordinate system.
%       h - handle to the figure
%
% Required Global Data Structures:
%       GLOBAL   - NRCIM
%       GLOBAL   - CST
%       GLOBAL   - MA
%
% Required Data Files:
%      none
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% static char rcsid[] = "$Id: nrCsPlotTriads.m,v 1.1 2012/07/24 18:40:27 roberts Exp $";
% INDENT-OFF*
% $Log: nrCsPlotTriads.m,v $
% Revision 1.1  2012/07/24 18:40:27  roberts
% Initial revision
%
% Revision 1.3  2005/07/29 22:20:45  roberts
% Now returns handle to figure
%
% Revision 1.2  2005/07/29 22:07:31  msmith
% Changed name of return parameter to be consistent with components.
%
% Revision 1.1  2004/10/26 19:39:58  roberts
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

function [xyuv h] = nrCsPlotTriads(surfnum)

global NRCIM;
global CST;
global MA;

numSegs = NRCIM.SegsPerSurf(surfnum);
x = zeros(numSegs,1);
y = x;
u = x;
v = x;

P1 = zeros(3,1);
P2 = zeros(3,1);
P3 = zeros(3,1);

% Create a list of the node numbers associated with the TRIADS in the
% surface
N=[CST.Meta([CST.Meta(:).NrcimSurf] == surfnum).NodeNums];
N = reshape(N',3,numSegs)';

for ii=1:numSegs
    % Get the node locations in the G coordinate system
    P1 = full(MA.NodePosG(N(ii,1),:));
    P2 = full(MA.NodePosG(N(ii,2),:));
    P3 = full(MA.NodePosG(N(ii,3),:));
    % transform the points to the Z coordinate system
    P1 = vlCsPMult(CST.ZG,P1');
    P2 = vlCsPMult(CST.ZG,P2');
    P3 = vlCsPMult(CST.ZG,P3');
    ave = (P1+P2+P3)/3;
    x(ii) = ave(1);
    y(ii) = ave(2);
    u(ii) = P3(1) - x(ii);
    v(ii) = P3(2) - y(ii);
end

figure;
h = quiver(x,y,u,v);
xyuv = [x;y;u;v];

% end of nrCsPlotTriads.m
        
        
    
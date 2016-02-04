% File: vlCsTriadTrans.m, NRCIM Toolbox
%
% Syntax: T = vlCsTriadTrans(P1,P2,P3)
%
% Discussion:
%   Given a set of three points P in global coordinates, 
%   the routine calculates the 4x4 transformation matrix T between 
%   the global coordinate system and this newly defined system.
%   The transform matrix is such that if P is a point in the local coordinate
%   system of the CS defined by (P1,P2,P3) then T * P is the point in
%   the global coordinate system.
%
%                    * P3                         Y Axis
%                   / \                           ^
%                  /   \                          |
%                 /  C  \ C = Centroid (Origin)   O - > X Axis
%                /       \
%           P1  *---------* P2 
%
%   The coordinate systems are constructed as follows:
%       The origin is at the centroid of the 3 points
%       The y axis points towards the upper point
%       The X axis points in the direction from P1 to P2
%       The Z axis forms a right handed coordinate system
%
% Input Parameters:
%   P1 - lower left point on triad (3 x 1 vector)
%   P2 - lower right point on triad (3 x 1 vector)
%   P3 - upper point on triad (3 x 1 vector)
%
%   P1,P2 and P3 are positions specified in the global coordinate system.
%
% Output Parameters:
%   T - 3x4 Transformation matrix (See Data Structures Document)
%
%   The transform is such that P_global = T * P_local, where P_local is
%   a point in the local coordinate system and P_global is the global
%   coordinate system (same CS as the input points P1,P2,P3 are defined in).
%
% Required Global Data Structures:
%       None
%
%
% Required Data Files:
%       None
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% static char rcsid[] = "$Id: vlCsTriadTrans.m,v 1.4 2004/12/06 23:14:42 msmith Exp $";
% INDENT-OFF*
% $Log: vlCsTriadTrans.m,v $
% Revision 1.4  2004/12/06 23:14:42  msmith
% Removed call to invert transformation matrix T.
% Updated comments to more explicitly define the meaning of matrix T.
%
% Revision 1.3  2004/12/01 06:27:13  roberts
% Inverted the transform
%
% Revision 1.2  2003/09/18 23:43:49  stretchn
% Updated header
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

function T=vlCsTriadTrans(P1,P2,P3)

if size(P1) ~= [3 1] | size(P2) ~= [3 1] | size(P3) ~= [3 1]
    error('Input points are not 3 x 1 vectors');
end

Z = vlCsUnit(cross(P2-P1,P3-P1));

C = (P1+P2+P3)/3;

Y = vlCsUnit(P3-C);

X = cross(Y,Z);

T = [X Y Z C];

% End of vlCsTriadTrans.m

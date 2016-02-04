% File: <vlOpMovetoSphere.m> , VLOT Toolbox
%
% Syntax: [RAYSOUT] = vlOpMovetoSphere(RAYSIN,d,c,zdir)
%
% Description:
%       This routine traces the rays along there paths in RAYSIN from 
%       there present location to the intersection with a sphere located 
%       on the Z axis at (x,y,z)=(0,0,d) with radius of curvature 
%       = 1/curvature. X,Y,Z and OPL of the ray are updated. The coordinate 
%       frame is unchanged.
%
%       The sign convention on curvature is the same as for Zemax:
%       Curvature is positive if the center of curvature is to the right
%       ( a positive distance along the z axis ) from the vertex of the sphere
%
%       All rays are tested for intersection with the sphere, if any ray
%       fails the routine is aborted with an error message.
%
%       See "Fundamental Optical Design", Kidger, 2002 P.51 for Algorithm
%
% Input Parameters:
%       RAYSIN    - (nx10 array ) Rays formatted in a (n,10) array format where 
%                   each row is ( s w i x y z l m n opl ). The number of rows
%                   is the number of rays.
%
%       d         - (scalar) The Z location of the vertex of the sphere in the local 
%                   coordinate system. The location of the vertex is
%                   (x,y,z)=(0,0,d). d is a physical distance.
%
%       c         - (scalar) The curvature of the reference sphere. curvature=1/Radius
%                   The sign convention on curvature is the same as for Zemax:
%                   Curvature is positive if the center of curvature is to the right
%                   ( a positive distance along the z axis ) from the vertex of the sphere
%
%       zdir      - (scalar, +1 or -1 ) The Z direction at the surface which 
%                   the rays are being traced to a sphere. If Z dir is positive, 
%                   the light is traveling in the direction of the Z axis. If Z
%                   is negitive, Z is traveling opposite the direction of the Z axis.
%
% Output Parameters:
%       RAYSOUT   - (nx10 array ) The rays intersected with the sphere. 
%                   X,Y,Z, and OPL have been updated in the calculation. 
%                   All other values are unchanged (S W I L M N ).
%
% Required Global Data Structures:
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
%   Copyright (c) 2003, John Pazder, National Research Council
%       Rev 1, Initial Release, May 27,2003  
%
% static char rcsid[] = "$Id: vlOpMovetoSphere.m,v 1.3 2004/11/10 22:22:23 pazderj Exp $";
% INDENT-OFF*
% $Log: vlOpMovetoSphere.m,v $
% Revision 1.3  2004/11/10 22:22:23  pazderj
% fix but with adding L to OPD now have zdir arg see bugaware 55 .
%
% Revision 1.2  2003/09/19 19:04:50  stukasa
% Header update.
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


function RAYSOUT = vlOpMovetoSphere(RAYSIN,d,c,zdir)

% check number of arguments
if  nargin ~= 4,  error('Wrong number of parameters, must enter RAYSIN,d and c ') ; end

% equate raysout with rays in to copy s w i L M N D and store Z of starting
% reference frame
RAYSOUT = RAYSIN ; 

% RAYS(:,n)   1 2 3 4 5 6 7 8 9 10
%             s w i x y z l m n D

% transform input rays to coordinate frame of vertex of sphere
RAYSIN(:,6) = RAYSIN(:,6)-d;

% compute intermediate values G and F for each ray ( G and F are vectors )
% G = N-c*(XL+YM+ZN)
G=RAYSIN(:,9)-c*(RAYSIN(:,4).*RAYSIN(:,7)+RAYSIN(:,5).*RAYSIN(:,8)+RAYSIN(:,6).*RAYSIN(:,9));
% F = c(X^2+Y^2+Z^2)-2*Z
F= c*(RAYSIN(:,4).^2+RAYSIN(:,5).^2+RAYSIN(:,6).^2)-2*RAYSIN(:,6);
% Q = G^2-cF 
Q=G.^2-c*F;
% test for beam missing the sphere i.e. Q<0
if any(Q<0), error('At least One Ray has missed the Sphere, check input parameters, aborting');end
% now compute length of ray from starting point to intersection of sphere
% for each ray
% P=F/(G+sqrt(G^2-cF))=F/(G+sqrt(Q))
P=F./(G+sqrt(Q));

% compute new intersections
% use Z in raysout so we are back in the original reference frame
% X'=X+LP
% Y'=Y+MP
% Z'=Z+NP
% OPL'=OPL+P
RAYSOUT(:,4)=RAYSOUT(:,4)+RAYSOUT(:,7).*P;
RAYSOUT(:,5)=RAYSOUT(:,5)+RAYSOUT(:,8).*P;
RAYSOUT(:,6)=RAYSOUT(:,6)+RAYSOUT(:,9).*P;
RAYSOUT(:,10)=RAYSOUT(:,10)+zdir*P;

% end of function
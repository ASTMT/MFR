% File: vlCsInvEulerXYZ.m, NRCIM Toolbox
% 
% Syntax: CST = vlCsInvEulerXYZ(Set)
%   
% Discussion:
%       This routine creates a 3x4 transformation matrix that corresponds
%       to a Euler XYZ transformation.
%       
%       The Euler X,Y,Z angle set is 
%	        (where c=cos, s=sin, A=alpha, B=beta, G=Gamma)
%                       |   cBcG            -cBsG           sB      |
%           Rxyz =      |   sAsBcG+cAsG     -sAsBsG+cAcG    -sAcB   |
%                       |   -cAsBcG+sAsG    cAsBsG+sAcG     cAcB    |
%
%       Note that a 4th row must be added before transformations are made.
%       This row is.
%
%                       |   0       0               0               1   |
%
%       Which is equivalent to
%           Rxyz = cstmult([trans([X Y Z]'),rot_x(A),rot_y(B),rot_z(G))
%
% Input Parameters:
%       Set - Vector of [X Y Z RX RY RZ] (R's in degrees)
%      
% Output Parameters:
%        CST - 4x4 transform matrix
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
% static char rcsid[] = "$Id: vlCsInvEulerXYZ.m,v 1.2 2003/09/18 23:43:49 stretchn Exp $";
% INDENT-OFF*
% $Log: vlCsInvEulerXYZ.m,v $
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

function CST = vlCsInvEulerXYZ(Set)

deg = 180/pi;

if size(Set,1) ~= 1 || size(Set,2) ~= 6
    error('Input must be a 1x6 array');
end

CST = zeros(3,4);
CST([1:3],4) = Set([1:3])';

% Convert angles into radians
Set([4:6]) = Set([4:6])/deg;

cA = cos(Set(4));   % alpha is the X rotation
cB = cos(Set(5));   % beta is the Y rotation
cG = cos(Set(6));   % gamma is the Z rotation
sA = sin(Set(4));
sB = sin(Set(5));
sG = sin(Set(6));

CST(1,[1:3]) = [cB*cG               -cB*sG              sB      ];
CST(2,[1:3]) = [sA*sB*cG+cA*sG      -sA*sB*sG+cA*cG     -sA*cB  ];
CST(3,[1:3]) = [-cA*sB*cG+sA*sG     cA*sB*sG+sA*cG      cA*cB   ];

% End of vlCsInvEulerXYZ
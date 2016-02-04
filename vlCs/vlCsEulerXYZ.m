% File: vlCsEulerXYZ.m, NRCIM Toolbox
%
% Syntax: Set = vlCsEulerXYZ(T)
%
% Discussion:
%       this routine finds the Euler X,Y,Z angle rotations required to
%       recreate the particular transformation
%       The Euler Z,Y,X angle set is 
%	        (where c=cos, s=sin, A=alpha, B=beta, G=Gamma)
%                       |   cBcG            -cBsG           sB      |
%           Rxyz =      |   sAsBcG+cAsG     -sAsBsG+cAcG    -sAcB   |
%                       |   -cAsBcG+sAsG    cAsBsG+sAcG     cAcB    |
%
%   To use the results, the correct order of transformation is ...
%       First, translate the original coordinate system by X,Y,Z
%       Then, rotate around the transformed X axis by RX degrees
%       Then, rotate about the transformed Y axis by RY degrees
%       Then, rotate about the transformed Z axis by RZ degrees
%           (I.E. these are Euler X,Y,Z rotations)
% Input Parameters:
%       T - 3x4 transform matrix (See CST in Data Structures Documentation)
% Output Parameters:
%       Set - Vector of [X Y Z RX RY RZ] (R's in degrees)
%
% Required Global Data Structures:
%       None
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
% static char rcsid[] = "$Id: vlCsEulerXYZ.m,v 1.2 2003/09/18 23:43:49 stretchn Exp $";
% INDENT-OFF*
% $Log: vlCsEulerXYZ.m,v $
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


function Set = vlCsEulerXYZ(T)

if size(T,1) ~= 3 || size(T,2) ~= 4
    error('Input must be a 3x4 array');
end

cB = sqrt(T(2,3)^2+T(3,3)^2);
B = atan2(T(1,3),cB);
if cB < eps
    A = 0;  % alpha can be chosen as 0
    G = atan2(T(2,1),T(2,2));
else
    G = atan2(-T(1,2),T(1,1));
    A = atan2(-T(2,3),T(3,3));
end
deg = 180/pi;
Set = [T(1,4) T(2,4) T(3,4) A*deg B*deg G*deg];

% End of vlCsEulerXYZ
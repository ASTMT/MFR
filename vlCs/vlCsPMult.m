%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: vlCsPMult.m, NRCIM Toolbox
%
% Syntax: pout = vlCsPMult(cst,pin)
%
% Discussion:
%   This transforms a point to a new coordinate system
%
% Input Parameters:
%       cst - CST transform matrix
%       pin - point [X Y Z]'
%
% Output Parameters:
%       pout - transformed point
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
% static char rcsid[] = "$Id: vlCsPMult.m,v 1.3 2004/10/01 21:27:31 msmith Exp $";
% INDENT-OFF*
% $Log: vlCsPMult.m,v $
% Revision 1.3  2004/10/01 21:27:31  msmith
% Removed unnecessary addition of fourth row of matrix (also eliminates
% the need to remove the fourth element of the multiplication result).
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

function pout = vlCsPMult(cst,pin)

if size(cst) ~= [3 4]
    error('Input CST matrix must have size 3 x 4');
end

if size(pin) ~= [3 1]
    error('Input point must have size 3 x 1');
end

% Since the cst matrix is being used to transform a point in 3-space to a new
% point in 3-space, we do not add the [ 0 0 0 1 ] row to the cst matrix but we
% do expand the pin vector by adding 1 (so that it is consistent with the 
% multiplication by the 3 x 4 matrix).

pout = cst * [ pin ; 1];

% End of vlCsPMult.m

% File: vlCsMult.m, NRCIM Toolbox
%
% Syntax: CST = vlCsMult(varargin)
%
% Discussion:
%   Multiplies several 3x4 CST coordinate transformations together. Takes
%   care of adding and removing the fourth row [0 0 0 1] of the
%   transformation matrix.
%
% Input Parameters:
%       List of CST's in the order they are to be multiplied
%
% Output Parameters:
%       CST - 3x4 transform matrix
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
% static char rcsid[] = "$Id: vlCsMult.m,v 1.2 2003/09/18 23:43:49 stretchn Exp $";
% INDENT-OFF*
% $Log: vlCsMult.m,v $
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

function CST = vlCsMult(varargin)
   
% get # of input arguments
N = nargin;

CST = eye(4);

for(i=1:N)
    CST = CST*[varargin{i}; 0 0 0 1];
end

CST = CST(1:3,1:4);

% End of vlCsMult.m
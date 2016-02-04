% File: vlCsTrans.m, NRCIM Toolbox
%
% Syntax: CST = vlCsTrans(P)
%
% Discussion:
%   This routine generates a CST data structure representing a translation
%   by the offset P.  CST is a 3x4 matrix as described in the Data
%   Structures document.
%
% Input Parameters:
%       P - [X Y Z]' offset
%
% Output Parameters:
%       CST - transform matrix
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
% static char rcsid[] = "$Id: vlCsTrans.m,v 1.2 2003/09/18 23:43:49 stretchn Exp $";
% INDENT-OFF*
% $Log: vlCsTrans.m,v $
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

function CST = vlCsTrans(P)

CST = [eye(3,3) P];
 
% End of vlCsTrans.m



% File: <vlLmRunLom.m>
%
% Syntax: [ W ] = vlLmRunLom(SP)
%
% Description:
%       Runs the Linear Optics Model
%
% Input Parameters:
%       SP - the perturbations 
%
% Output Parameters:
%       W - the wavefront map, N^2 x 1 vector
%
% Required Global Data Structures:
%       OC
%       IM
%
% Required Data Files:
%              

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% static char rcsid[] = "$Id: vlLmRunLom.m,v 1.4 2003/11/24 19:50:18 stretchn Exp $";
% INDENT-OFF*
% $Log: vlLmRunLom.m,v $
% Revision 1.4  2003/11/24 19:50:18  stretchn
% Now outputs dummy OPD for Zemax case, so switch inputs will be of the same
% length
%
% Revision 1.3  2003/06/26 22:28:45  mckenzie
% updated headers and comments
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

function OPD = vlLmRunLom(SP)

global IM
global OC

if strcmp(IM.OpticalEngine,'LOM')
    OPD = OC.LOM.OPDzero + OC.LOM.dWdP*SP;
%These statements are necessary so that input vectors to the optical engine
%switch will always be the same length
else
    global RAYSIN;
    OPD = zeros(length(RAYSIN),1);
end
% end <vlLmRunLom.m>

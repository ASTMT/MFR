% File: <nrCimCP.m>
%
% Syntax: CP = nrCimCP(CstIdx,Point)
%
% Description:
%       Returns the position of a triad point in its local C coordinate.
%
% Input Parameters:
%       CstIdx      - (scalar integer) Surface number as per the order in
%                       the CST.Meta data structure.
%       Point       - (scalar integer) Point number in triad
%
% Output Parameters:
%       CP        - (3x1) Position of point in its local C coordinate system
%
% Required Global Data Structures:
%       MA, CST
%
% Required Data Files:
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: nrCimCP.m,v 1.1 2008/04/30 19:35:57 roberts Exp $
%
% INDENT-OFF*
% $Log: nrCimCP.m,v $
% Revision 1.1  2008/04/30 19:35:57  roberts
% Initial revision
%
% INDENT-ON*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%           Herzberg Institute of Astrophysics                  %%%%%
%%%%%%      Astronomy Technology Research Group - Victoria           %%%%%
%
% (c) <2008>				        (c) <2008>
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

function CP = nrCimCP(CstIdx,Point)

global MA
global CST

if Point < 1 || Point > 3
    error('Point exceeds bounds');
end
    
% get the index into the node list
NdIdx = CST.Meta(CstIdx).NdLstStart + Point - 1;

% get the position of the node in the G coordinate system
GP = MA.NodeLst(NdIdx,2:4)';

% transform to the C coordinate position

CP = vlCsPMult(vlCsMult(vlCsInv(CST.ZC(:,:,CstIdx)),CST.ZG),GP);

% End of nrCimCP
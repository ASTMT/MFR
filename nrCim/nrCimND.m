% File: <nrCimND.m>
%
% Syntax: ND = nrCimND(CstIdx,Point,Y)
%
% Description:
%       Provides the nodal displacement in the nodal coordinate system for
%       a given surface, triad point and Y vector.
%
% Input Parameters:
%       CstIdx    - (scalar integer) CST index number
%       Point     - (scalar integer) Point number.  Must be 1 for '1Node'
%                   and 1 to 3 for 'Triad' surfaces.
%
% Output Parameters:
%       ND       - (3 x 1) Nodal displacement vector (Txyz) for 'Triad' type points
%                - (6 x 1) Nodal displacement vector (Txyz,Rxyz) for '1Node' Type Points
%
% Required Global Data Structures:
%       CST
%
% Required Data Files:
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: nrCimND.m,v 1.2 2009/06/25 22:12:29 roberts Exp $
%
% INDENT-OFF*
% $Log: nrCimND.m,v $
% Revision 1.2  2009/06/25 22:12:29  roberts
% fixed error in indexing into Y - now uses YMAP
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

function ND = nrCimND(CstIdx,Point,Y)

global CST

% Get the starting row number in the Y vector

%fprintf('CstIdx,Point = (%d,%d)\n',CstIdx,Point);

switch CST.Meta(CstIdx).AType
    case 'Triad'
        if Point > 3 || Point < 1
            error('Point must be between 1 and 3 for Triad CSTs');
        end
        rowstart = 1 + (Point-1)*3;
        ND = Y(CST.Meta(CstIdx).YMap(rowstart:rowstart+2));

    case '1Node'
        if Point > 1
            error('Point must equal 1 for 1Node CSTs');
        end
        ND = Y(CST.Meta(CstIdx).YMap);
    otherwise
        error('Bad AType')
end

% End of nrCimND
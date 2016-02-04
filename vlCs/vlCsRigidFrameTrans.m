% File: vlCsRigidFrameTrans.m, NRCIM Toolbox
%
% Syntax: [CSTBBP,EBBP] = vlCsRigidFrameTrans(CSTZA,CSTZAP,CSTZB)
%
% Discussion:
%   Given 2 coordinate systems, A and B, fixed to a rigid body, consider the case of 
%   the rigid body moving to a position described by a new CS A'.  What is
%   the corresponding CS B',and what are the Euler angles to get from B to B'?
%
% Input Parameters:
%       CSTZA - Transform between fixed CS Z and A
%       CSTZAP - Transform between fixed CS Z and A'
%       CSTZB - Transform between fixed CS Z and B
%       
% Output Parameters:
%       CSTBBP - Transform between B and B'
%       EXYZBBP - Euler X,Y,Z,Rx,Ry,Rz angles between B and B'
%
% Required Global Data Structures:
%       None
%
% Required Data Files:
%       None
%       

%
% Verification
%
% CSTZA = vlCsTrans([0 10 0]');
% CSTZAP = vlCsMult(vlCsTrans([0 11 0]'),vlCsRotZ(90));
% CSTZB = vlCsMult(vlCsTrans([0 12 0]'),vlCsRotZ(180));
% 
% % by inspection EXYZBBP should be X=2,Y=1,Z=0,Rx=0,Ry=0,Rz=90
% 
% [CSTBBP,EBBP] = vlCsRigidFrameTrans(CSTZA,CSTZAP,CSTZB);

% Results
% CSTBBP =
%     0.0000   -1.0000         0    2.0000
%     1.0000    0.0000         0    1.0000
%          0         0    1.0000         0
% EBBP =
%     2.0000    1.0000         0         0         0   90.0000
%

%
% Revision History
%
% static char rcsid[] = "";
% INDENT-OFF*
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
function [CSTBBP,EXYZBBP] = vlCsRigidFrameTrans(CSTZA,CSTZAP,CSTZB)

CSTAB = vlCsMult(vlCsInv(CSTZA),CSTZB);
CSTAPBP = CSTAB;    % by definition since they are on a rigid body
CSTBBP = vlCsMult(vlCsInv(CSTAB),vlCsInv(CSTZA),CSTZAP,CSTAPBP);
EXYZBBP = vlCsEulerXYZ(CSTBBP);


% end of vlCsRigidFrameTrans





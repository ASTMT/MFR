% File: <vlSciActStroke.m>, NRCIM Toolbox
%
% Syntax: stroke = vlSciActStroke(NomP, NewP)
%
% Description:
%       This routine calculates the sag of a conic surface.
%   
% Input Parameters:
%       NomP - [1x3 vector] nominal actuator position in Z CS
%       NewP - [1x3 vector] new actuator position in Z CS
%
% Output Parameters:
%       stroke - Required actuator stroke for the change of actuator
%       position from Nom to New in C CS (aligned with actuator direction)
%
% Required Global Data Structures:
%       none
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
%   Copyright (c) 2010, Scott Roberts, National Research Council 
%
% static char rcsid[] = "$Id: vlSciActStroke.m,v 1.3 2012/10/24 23:23:10 roberts Exp $";
% INDENT-OFF*

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

function stroke = vlSciActStroke(NomP, NewP)

% Calculates actuator stroke based on a nominal position and a new position
global OC;
global CST;
global TMTCS;

% Get the radius in the Z CS plane for each input point
NomR = norm(NomP(1:2));
NewR = norm(NewP(1:2));

% Calculate the sag and slope at each position
[NomZ NomSlope] = vlOpSag(1/OC.PrimRadCurv,OC.PrimCC,NomR);
[NewZ NewSlope] = vlOpSag(1/OC.PrimRadCurv,OC.PrimCC,NewR);

% Get the nominal slope.  Z stroke increases by 1/cos(theta)
theta = atan(NomSlope);

% Actuator stroke is in C coordinates, positive in C Z-direction
% Sag is positive for increasing radii
% Need to adjust the sag in the direction of actuator motion
% Sag is positive in TMT M1CRS direction
% Check Z CS Z axis direction

% Are Z and C aligned?
ZCsign = 1;
if CST.ZC(3,3,1) < 0
    ZCsign = -1;
end

% Are the Sag direction and C aligned?
if abs(TMTCS.ZM1(3,3)) ~= 1
    error('Z axis of Z and M1CRS axes are not aligned');
end
SagCsign = TMTCS.ZM1(3,3)*ZCsign;

% The stroke is the difference between the change in sag and the exi

% The sag is required stroke so should be added
strokeZSag = SagCsign*(NewZ-NomZ);
% strokeZP is an error so needed stroke is -strokeZP
strokeZP = ZCsign*(NewP(3)-NomP(3));

% strokeZ
strokeZ = strokeZSag - strokeZP;

stroke = strokeZ * cos(theta);
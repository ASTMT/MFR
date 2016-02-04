% File: Set = vlOpDecComp(TZI,Dx,Dy,PivotDist)
%
% Calculates the required motions to best compensate a segment decenter 
% in the C coordinate system X,Y direction. Also provides the delta X and 
% Y displacement due to the coupled motion of tilt and decenter.
% Input Parameters:
%       TZI         - (scalar) the ZMNonSeq coordinate transform to the seg
%       exyz        - (1x6 vector) EulerXYZ (dx,dy,dz,rx,ry,rz) motions of
%                   the segment that describe the M to M' transform.
%
% Output Parameters
%       Set     - [6 x 1]. The full EulerXYZ angle set in the form
%                   dx,dy,dz,rx,ry,rz

% The SSA design has a "pivot point" below the segment that results in a
% decenter of the segment for a tilt. The pivot point location is given by
% OC.PivotDist.  The Dx and Dy that are calculated from the Y vector are
% due to mirror cell deformations and are not because of the pivot point
% location.  However, when correcting the segment tip/tilt due to mirror
% cell location there is a motion in Dx and Dy. The adjustment to tilt the
% mirror to its nominally aligned position is calculated. There is an
% additional motion required for the offset in sag, such that a small tilt
% of the segment is introduced for a new off-axis position and this incurs
% an additional decenter. Some iteration would be required to get the exact
% solution, but the result here is sufficient for the purposes of the MFR.
%
% If the segment has rotated by Rx,Ry and decentered by Dx,Dy then the
% adjustment should be as follows (where pivotdist is a positive value for
% below segment i.e. PivotDist is a positive number)
%
% DeltaDx = -sin(Ry)*pivotdist DeltaDy = sin(Rx)*pivotdist

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% static char rcsid[] = "$Id: vlOpDecComp.m,v 1.5 2012/10/26 20:38:09 roberts Exp roberts $";
% INDENT-OFF*
% $Log: vlOpDecComp.m,v $
% Revision 1.5  2012/10/26 20:38:09  roberts
% EXYZ return values are same as input except for Dx and Dy.  Updated comments.
%
% Revision 1.4  2011/03/02 16:49:09  vlotim
% Commented out fprintf statements - Scott
%
% Revision 1.3  2010/05/19 20:23:51  roberts
% updated for pivot point on segment and related decenter for tilt of segment
%
% Revision 1.2  2004/11/30 04:49:12  roberts
% now returns the full Euler angle set as the 4th return value
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

function Set = vlOpDecComp(TZI,exyz)

global OC;

% check number of arguments
if  nargin ~= 2,  error('Wrong number of parameters, must enter TZI,Set') ; end

PivotDist = OC.PivotDist;

% PivotDist   - (scalar) distance of segment pivot point below M1 segment
% surface. This distance is in the positive Z direction in the segment M
% coordinate system since M points away from the sky, and the pivot point
% is below the mirror surface.

Dx = exyz(1);
Dy = exyz(2);
Rx = exyz(4);
Ry = exyz(5);

%fprintf('\nX,Y = %.2f,%.2f\t',TZI(1,4),TZI(2,4));
%fprintf('Idxy = %.2e,%.2e\t',Dx,Dy);

% Adjust Dx and Dy for compensating the Rx and Ry of the segment motion
Dx = Dx - sin(Ry*pi()/180)*PivotDist;
Dy = Dy + sin(Rx*pi()/180)*PivotDist;
%fprintf('Adxy = %.2e,%.2e\t',Dx,Dy);

% calculate the displaced position of the segment in Z coordinates
PD = vlCsPMult(TZI,[Dx Dy 0]');
SegXD = PD(1);
SegYD = PD(2);

% convert to polar 
[thetaD,rhoD] = cart2pol(SegXD,SegYD);
[thetaI,rhoI] = cart2pol(TZI(1,4),TZI(2,4));

Rz = (thetaD - thetaI)*180/pi;
%fprintf('Rz = %.2e\n',Rz);

% want orientation of Y axis, so subtract pi/2
thetaD = thetaD - pi/2;

% compute sag and radial slope at the displaced position
[zD, dzdrD] = vlOpSag(1/OC.PrimRadCurv,OC.PrimCC,rhoD);

% calculate transform from the new position back to Z
TZD = vlCsMult(vlCsTrans([SegXD SegYD -zD]'),vlCsRotZ(thetaD*180/pi),vlCsRotX(-atan(dzdrD)*180/pi));

% calculate the transform from TZI to TZD. It is okay that the z rotation
% of the two CS don't match, since we only want Rx,Ry values and we're
% calclulating Euler XYZ angles.
TID = vlCsMult(vlCsInv(TZI),TZD);

% calculate the rotation angles
Set = vlCsEulerXYZ(TID);
Rx = Set(4);
Ry = Set(5);

% Calculate the new decenter after the tip-tilt motion (because of pivot distance)
Set(1) = Set(1) - sin(Ry*pi()/180)*PivotDist;
Set(2) = Set(2) + sin(Rx*pi()/180)*PivotDist;

% Z piston correction is by definition zero
Set(3) = 0;

% Calculated Rx and Ry values are necessary to keep in sp so that OPD plots
% compensate tilts with decenters.

% The two coordinate systems are rotated differently around Z so replace
% the Rz from the set with the value calculated above which is the clocking
% of the segment
%Set(6) = Rz;
Set(6) = exyz(6);


%fprintf('Fdxy = %.2e,%.2e\n',Set(1),Set(2));
% end of vlOpDecComp





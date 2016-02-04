% File: <nrSciCalcMotion2Cases.m>
%
% Syntax: motions = nrSciCalcMotion2Cases(basisRES,pertRES,ElevAng,fid)
%
% Description:
% This function reports the motion of the M1, M2, M3 and Focal Plane
% between two MFR outputs.  The first input argument,Basis, is the case
% where all optics are assumed to be in alignment.  The second case,
% Perturbed, is the case where motions relative to the Basis case are
% reported.
%
% Input Parameters:
%       basisRESFILE - MFR output RES data structure
%       pertRESFILE - MFR output RES data structure
%       ElevAng - Elevation angle of pert case minus the basis case (in
%       degrees)
%       fid - file ID to write to
%
% Output Parameters:
%       motions - 
%
% Required Global Data Structures:
%       GLOBAL1   - (Optional Comments, listing of fields)
%       GLOBAL2   - (Optional Comments, listing of fields)
%
% Required Data Files:
%       File1.xxx - Description
%       File2.xxx - Description
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: nrSciCalcMotion2Cases.m,v 1.11 2012/10/29 22:53:00 roberts Exp $
%
% INDENT-OFF*
% $Log: nrSciCalcMotion2Cases.m,v $
% Revision 1.11  2012/10/29 22:53:00  roberts
% tweaks to report format
%
% Revision 1.10  2012/10/29 22:48:17  roberts
% REPORT now says EulerXYZ angles
%
% Revision 1.9  2012/10/29 22:19:44  roberts
% updated report language to clarify CS and order of Euler angles and translations
%
% Revision 1.8  2012/10/15 15:52:12  roberts
% updated FOS to M3 transformation calculation and reporting
%
% Revision 1.7  2012/10/09 22:25:54  roberts
% updated outputs and formatting of report
%
% Revision 1.6  2012/10/05 21:07:50  roberts
% fixed report formatting
%
% Revision 1.5  2012/10/02 22:09:23  roberts
% Corrected the reported names of the Basis and Perturbed cases.
%
% Revision 1.4  2012/07/24 20:00:32  roberts
% major update to allow calling from vlSciWriteResults.m
%
% Revision 1.3  2011/02/18 00:48:01  roberts
% Added additional results
%
% Revision 1.1  2010/02/22 01:43:06  roberts
% Initial revision
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

function GRAVRES = nrSciCalcMotion2Cases(basisRES,pertRES,ElevAng,fid)

% Coordinate systems are:
%   Aligned relative to M1 actuator fit
%   Original relative to M1 nominal
%   Displaced at uncorrected position due to FEA load

% Input cases are:
%   Basis: the case at which all optics are perfectly aligned
%   Pert: the perturbed case to be compared with basis

%basis = load(basisRESFILE);
%delta = load(pertRESFILE);

basis = basisRES;
delta = pertRES;




%% Nominal Transforms
% Calculate nominal TMT CS transforms
CST_Z_TMTM1 = basis.CST_Z_TMTM1;
CST_Z_TMTM2 = basis.CST_Z_TMTM2;
CST_Z_TMTM3 = basis.CST_Z_TMTM3;
CST_Z_TMTFOS = basis.CST_Z_TMTFOS;

T_M1O_M2O = vlCsMult(vlCsInv(CST_Z_TMTM1),CST_Z_TMTM2);
T_M1A_M2A = T_M1O_M2O;
T_M1O_M3O = vlCsMult(vlCsInv(CST_Z_TMTM1),CST_Z_TMTM3);
T_M1A_M3A = T_M1O_M3O;
T_M1O_FOSO = vlCsMult(vlCsInv(CST_Z_TMTM1),CST_Z_TMTFOS); % FOS coordinate system rotates about Z as the telescope rotates
T_M1A_FOSA = T_M1O_FOSO;
T_FOSA_M3A = vlCsMult(vlCsInv(T_M1A_FOSA),T_M1A_M3A);

GRAVRES.T_M1A_M2A = T_M1A_M2A;
GRAVRES.T_M1O_M2O = T_M1O_M2O;
GRAVRES.T_M1A_M3A = T_M1A_M3A;
GRAVRES.T_M1O_M3O = T_M1O_M3O;
GRAVRES.T_M1A_FOSA = T_M1A_FOSA;
GRAVRES.T_M1O_FOSO = T_M1O_FOSO;

%% M1 Primary Mirror

% Get the transform between M1O and M1A for the basis and perturbed cases

T_M1OB_M1AB = vlCsInvEulerXYZ(basis.EXYZM1P);
T_M1OP_M1AP = vlCsInvEulerXYZ(delta.EXYZM1P);

% Note that M1OB and M1OP vary by the difference in elevaton angle,
% therefore:
T_M1OB_M1OP = vlCsRotX(ElevAng);

% Now we can calculate how much the aligned M1's have moved between cases
T_M1OB_M1AP = vlCsMult(T_M1OB_M1OP,T_M1OP_M1AP);
T_M1AB_M1AP = vlCsMult(vlCsInv(T_M1OB_M1AB),T_M1OB_M1OP,T_M1OP_M1AP);
EXYZ_M1AB_M1AP = vlCsEulerXYZ(T_M1AB_M1AP);

GRAVRES.T_M1OB_M1AP = T_M1OB_M1AP;
GRAVRES.T_M1OB_M1AB = T_M1OB_M1AB;
GRAVRES.T_M1OP_M1AP = T_M1OP_M1AP;
GRAVRES.T_M1AB_M1AP = T_M1AB_M1AP;
GRAVRES.EXYZ_M1AB_M1AP = EXYZ_M1AB_M1AP;

%% M2

% Get the transform between M2A and M2D for the basis and perturbed cases
T_M2AB_M2DB = vlCsInvEulerXYZ(basis.EULM2A_M2D);
T_M2AP_M2DP = vlCsInvEulerXYZ(delta.EULM2A_M2D);

% Now calculate the transform between M2D basis and M2D perturbed
% This is how much the M2 has moved between the two cases
% Use the fact that M2AB will be the same position as M2AP
T_M2AB_M2DP = T_M2AP_M2DP;
T_M2DB_M2DP = vlCsMult(vlCsInv(T_M2AB_M2DB),T_M2AB_M2DP);

% Find the displacement of the M2 origin in M1 coordinates
% Assume that the total displacement occurs from a perfectly aligned base
% case.  I.E. Transform T_M2DB_M2DP into M1A coordinates
T_M2AP_M2DTotal = T_M2DB_M2DP;
T_M1AP_M2DP = vlCsMult(T_M1A_M2A,T_M2AP_M2DTotal);

% Calculate delta motion of M2 origin in M1 coordinates
D_XYZ_M2D_M1A = T_M1AP_M2DP(1:3,4)' - T_M1A_M2A(1:3,4)';

% Write results
GRAVRES.T_M2AB_M2DB = T_M2AB_M2DB;
GRAVRES.T_M2AP_M2DP = T_M2AP_M2DP;
GRAVRES.T_M2DB_M2DP = T_M2DB_M2DP;
GRAVRES.T_M2AP_M2DP = T_M2AP_M2DP;
GRAVRES.T_M1AP_M2DP = T_M1AP_M2DP;

GRAVRES.D_XYZ_M2D_M1A = D_XYZ_M2D_M1A;

EXYZ_M2DB_M2DP = vlCsEulerXYZ(T_M2DB_M2DP);
EXYZ_M1AP_M2DP = vlCsEulerXYZ(T_M1AP_M2DP);
GRAVRES.EXYZ_M2DB_M2DP = EXYZ_M2DB_M2DP;
GRAVRES.EXYZ_M1AP_M2DP = EXYZ_M1AP_M2DP;

%% M3
% Same as M2 for M3
T_M3AB_M3DB = vlCsInvEulerXYZ(basis.EULM3A_M3D);
T_M3AP_M3DP = vlCsInvEulerXYZ(delta.EULM3A_M3D);

% The M3 CS rotates as the telescope tilts in elevation
T_M3AB_M3AP = vlCsMult(vlCsRotX(45),vlCsRotZ(-ElevAng),vlCsRotX(-45));

T_M3AB_M3DP = vlCsMult(T_M3AB_M3AP,T_M3AP_M3DP);
T_M3DB_M3DP = vlCsMult(vlCsInv(T_M3AB_M3DB),T_M3AB_M3DP);

% Calculate M3 motion relative to FOS Basis case
T_FOSB_M3DP = vlCsMult(T_FOSA_M3A,T_M3DB_M3DP);

T_M3AP_M3DTotal = T_M3DB_M3DP;
T_M1AP_M3DP = vlCsMult(T_M1A_M3A,T_M3AP_M3DTotal);

%D_XYZ_M3D_M1A = T_M1AP_M3DP(1:3,4)' - T_M1A_M3A(1:3,4)';

EXYZ_M3DB_M3DP = vlCsEulerXYZ(T_M3DB_M3DP);
EXYZ_M1AP_M3DP = vlCsEulerXYZ(T_M1AP_M3DP);

GRAVRES.T_M3DB_M3DP = T_M3DB_M3DP;
GRAVRES.EXYZ_M3DB_M3DP = EXYZ_M3DB_M3DP;
GRAVRES.T_M3AB_M3DB = T_M3AB_M3DB;
GRAVRES.T_M3AP_M3DP = T_M3AP_M3DP;
%GRAVRES.D_XYZ_M3D_M1A = D_XYZ_M3D_M1A;
GRAVRES.T_M3AP_M3DP = T_M3AP_M3DP;
GRAVRES.T_M1AP_M3DP = T_M1AP_M3DP;
GRAVRES.EXYZ_M1AP_M3DP = EXYZ_M1AP_M3DP;
GRAVRES.T_FOSB_M3DP = T_FOSB_M3DP;
GRAVRES.EXYZ_FOSB_M3DP = vlCsEulerXYZ(T_FOSB_M3DP);

%% FOS
% Same for FOS
T_FOSAB_FOSDB = vlCsInvEulerXYZ(basis.EULFOSA_FOSD);
T_FOSAP_FOSDP = vlCsInvEulerXYZ(delta.EULFOSA_FOSD);

% Have to include the elevation angle change in this case since we are
% rotating the coordinate system in relation to to the structure that is
% deflecting.

T_FOSAB_FOSAP = vlCsRotZ(ElevAng);
T_FOSAB_FOSDP = vlCsMult(T_FOSAB_FOSAP,T_FOSAP_FOSDP);
T_FOSDB_FOSDP = vlCsMult(vlCsInv(T_FOSAB_FOSDB),T_FOSAB_FOSDP);
T_FOSDB_FOSDP_NoElev = vlCsMult(T_FOSDB_FOSDP,vlCsRotZ(-ElevAng));

%T_FOSAP_FOSDTotal = T_FOSDB_FOSDP;
%T_M1AP_FOSDP = vlCsMult(T_M1A_FOSA,T_FOSAP_FOSDTotal);
%D_XYZ_FOSD_M1A = T_M1AP_FOSDP(1:3,4)' - T_M1A_FOSA(1:3,4)';

EXYZ_FOSDB_FOSDP = vlCsEulerXYZ(T_FOSDB_FOSDP);
%EXYZ_M1AP_FOSDP = vlCsEulerXYZ(T_M1AP_FOSDP);

GRAVRES.T_FOSAB_FOSDB = T_FOSAB_FOSDB;
GRAVRES.T_FOSAP_FOSDP = T_FOSAP_FOSDP;
GRAVRES.T_FOSDB_FOSDP = T_FOSDB_FOSDP;
GRAVRES.EXYZ_FOSDB_FOSDP = EXYZ_FOSDB_FOSDP;
GRAVRES.T_FOSAP_FOSDP = T_FOSAP_FOSDP;
%GRAVRES.D_XYZ_FOSD_M1A = D_XYZ_FOSD_M1A;
GRAVRES.EULFOSO_FOSD = delta.EULFOSO_FOSD;

%% FOS to M3 Transformation
T_FOSDB_M3DP = vlCsMult(T_FOSA_M3A,T_M3DB_M3DP);
T_FOSDB_FOSDP_NoElev = vlCsMult(T_FOSDB_FOSDP,vlCsRotZ(-ElevAng));
T_FOSDP_NoElev_M3DP = vlCsMult(vlCsInv(T_FOSDB_FOSDP_NoElev),T_FOSDB_M3DP);

GRAVRES.T_FOSDP_NoElev_M3DP = T_FOSDP_NoElev_M3DP;
EXYZ_T_FOSDP_NoElev_M3DP = vlCsEulerXYZ(T_FOSDP_NoElev_M3DP);
GRAVRES.EXYZ_T_FOSDP_NoElev_M3DP = EXYZ_T_FOSDP_NoElev_M3DP;


%% FOS to M1 Transformations
%GRAVRES.T_M1AP_FOSDP = T_M1AP_FOSDP;
%GRAVRES.EXYZ_M1AP_FOSDP = EXYZ_M1AP_FOSDP;

%% Print Results


fprintf(fid,'\n\nMFR Output of Motions for Gravity Loading between two elevation angles\n\n');

fprintf(fid,'\tDefinitions:\n');
fprintf(fid,'\t\tBasis Case: the case at which all optics are perfectly aligned (usually at zenith)\n');
fprintf(fid,'\t\tPerturbed Case: the perturbed case to be compared with basis (usually off zenith elevation angle)\n\n');

fprintf(fid,'\tCoordinate Systems:\n');
fprintf(fid,'\t\tOriginal CS: as per OAD document\n');
fprintf(fid,'\t\tAligned CS: relative to M1 actuator fit (M1 Aligned)\n');
fprintf(fid,'\t\tDisplaced CS: position under displacement in FEA model\n\n');

fprintf(fid,'\tInput Cases:\n');
fprintf(fid,'\t\tBasis case is %s\n',basis.outFileRoot);
fprintf(fid,'\t\tPerturbed case is %s\n',delta.outFileRoot);
fprintf(fid,'\t\tThe difference in elevation angle between the two cases is %d degrees\n',ElevAng);

fprintf(fid,'\n\tActuator Stroke Calculation Information\n');
fprintf(fid,'\t\tMethod of Actuator Stroke Minimization [ActFit] is: %s\n',delta.ActFit);
if delta.NomAct
    fprintf(fid,'\t\tActuator Stroke is Calculated Relative to another load case\n');
    fprintf(fid,'\t\tBasis Load Case File [NomActPosFile] is : %s\n',delta.NomActPosFile);
else
    fprintf(fid,'\t\tActuator Stroke is Calculated Relative to no load case (zero thermal & gravity)\n');  
end
fprintf(fid,'\n');

fprintf(fid,'\tOutput displacements and angles are in meters and degrees\n\n');

% Coordinate System Output
fprintf(fid,'\tNominal Coordinate System Transformations\n');

fprintf(fid,'\t\tM1 [Aligned CS] to M2 [Aligned CS]\n');
fprintf(fid,'\t\t\tDisplacements, then EulerXYZ angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3f,%.3f,%.6f,%.3f,%.3f,%.3f\n',vlCsEulerXYZ(GRAVRES.T_M1O_M2O));
fprintf(fid,'\t\tM1 [Aligned CS] to M3 [Aligned CS]\n');
fprintf(fid,'\t\t\tDisplacements, then EulerXYZ angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n',vlCsEulerXYZ(GRAVRES.T_M1O_M3O));
fprintf(fid,'\t\tM1 [Aligned CS] to FOS [Aligned CS]\n');
fprintf(fid,'\t\t\tDisplacements, then EulerXYZ angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n',vlCsEulerXYZ(GRAVRES.T_M1O_FOSO));

% Primary Mirror Output (M1)
fprintf(fid,'\n\tPrimary Mirror (M1) Motions\n');
fprintf(fid,'\t\tTransform to M1 [Perturbed Case, Aligned CS] from M1 [Perturbed Case, Original Coordinates]\n');
fprintf(fid,'\t\t\tDisplacements, then EulerXYZ angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',GRAVRES.EXYZ_M1AB_M1AP);

% Secondary Mirror Output (M2)
fprintf(fid,'\n\tSecondary Mirror (M2) Motions\n');
fprintf(fid,'\t\tTransform to M2 [Perturbed Case, Displaced CS] from M2 [Perturbed Case, Aligned CS]\n');
fprintf(fid,'\t\t\tDisplacements, then EulerXYZ angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',GRAVRES.EXYZ_M2DB_M2DP);

fprintf(fid,'\t\tTransform to M2 [Perturbed Case, Displaced CS] from M1 [Perturbed Case, Aligned CS]\n');
fprintf(fid,'\t\t\tDisplacements, then EulerXYZ angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.6e,%.3e,%.3e,%.3e\n',GRAVRES.EXYZ_M1AP_M2DP);

fprintf(fid,'\t\tM2 [Perturbed Case, Displaced CS] displacement in M1 [Perturbed Case, Aligned CS]\n');
fprintf(fid,'\t\t\t(DX,DY,DZ) is...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e\n',GRAVRES.D_XYZ_M2D_M1A);

% Tertiary Mirror Output (M3)
fprintf(fid,'\n\tTertiary Mirror (M3) Motions\n');
fprintf(fid,'\t\tTransform to M3 [Perturbed Case, Displaced CS] from M3 [Basis Case, Aligned CS]\n');
fprintf(fid,'\t\t\tDisplacements, then EulerXYZ angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',GRAVRES.EXYZ_M3DB_M3DP);

fprintf(fid,'\t\tTransform to M3 [Perturbed Case, Displaced CS] from M1 [Basis Case, Aligned CS]\n');
fprintf(fid,'\t\t\tDisplacements, then EulerXYZ angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.6e,%.3e,%.3e,%.3e\n',GRAVRES.EXYZ_M1AP_M3DP);

fprintf(fid,'\t\tTransform to M3 [Perturbed Case, Displaced CS] from FOS [Basis Case, Aligned CS]\n');
fprintf(fid,'\t\t\tDisplacements, then EulerXYZ angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.6e,%.3e,%.3e,%.3e\n',GRAVRES.EXYZ_FOSB_M3DP);

% Focal Surface Output (FOS)
fprintf(fid,'\n\tInstrument Focal Surface (FOS) Motions\n');
fprintf(fid,'\t\tMFR FOS Coordinate System rotates with elevation angle\n');
fprintf(fid,'\t\tTransform to FOS [Perturbed Case, Displaced CS] from FOS [Basis Case, Aligned CS]\n');
fprintf(fid,'\t\t\tDisplacements, then EulerXYZ angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',GRAVRES.EXYZ_FOSDB_FOSDP);

% fprintf(fid,'\t\tFOS [Perturbed Case, Displaced CS] relative to M1 [Basis Case, Aligned CS]\n');
% fprintf(fid,'\t\t\tDisplacements, then EulerXYZ angles (DX,DY,DZ,RX,RY,RZ) are...\n');
% fprintf(fid,'\t\t\t\t%.3e,%.3e,%.6e,%.3e,%.3e,%.3e\n',GRAVRES.EXYZ_M1AP_FOSDP);

% M3 Relative to the Focal Surface

fprintf(fid,'\n\tM3 Position Relative to the Instrument Focal Surface (FOS)\n');
fprintf(fid,'\t\tTransform to M3 [Perturbed Case, Displaced CS] from FOS [Perturbed Case, Displaced CS, no elevation angle rotation of CS]\n');
fprintf(fid,'\t\t\tDisplacements, then EulerXYZ angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.6f,%.3e,%.3e,%.3e\n',GRAVRES.EXYZ_T_FOSDP_NoElev_M3DP);

dM3FOS = GRAVRES.EXYZ_T_FOSDP_NoElev_M3DP(1:3) - T_FOSA_M3A(:,4)';
fprintf(fid,'\t\tM3 to FOS, delta position from aligned coordinates\n');
fprintf(fid,'\t\t\t(DX,DY,DZ) is...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e\n',dM3FOS);
% End of nrSciCalcMotion2Cases


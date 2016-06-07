% File: alignment.m
%
% Syntax: [T_M1M2, T_M1M3, T_M1FOS] = alignment(zenith_angle, silent (optional))
%
% Discussion:
%   This routine generates coordinate transforms from M1 to M2, M1 to M3
%   and and M1 to FOS. These represent the gravity deformation between the
%   angle that the telescope is aligned at and the telescope zenith angle.
%
% Input Parameters:
%       zenith_angle - the zenith angle that deformations are reported at
%       silent - true or false to print information [default - false]
%       
% Input Files:
%       This function expects to find a file named XXalignYY_Data.mat where
%       XX is the zenith angle and YY is the angle the telescope is aligned
%       at.  This file contains the MFR results that are used to calculate
%       the transforms.
%
% Output Parameters:
%       T_M1M2 - transform matrix from M1 to the displaced M2
%       T_M1M3 - transform matrix from M1 to the displaced M3
%       T_M1FOS - transform matrix from M1 to the displaced Focal Surface

function [T_M1M2, T_M1M3, T_M1FOS] = alignment(zenith_angle, varargin)

if nargin < 2 
    silent = false;
else
    silent = varargin{1};
end

%% Setup
format short
matfile = sprintf('%dalign30_Data.mat',zenith_angle);

if ~silent 
    fprintf('...Loading %s\n', matfile)
end

load(matfile)
align_ang = sscanf(RES.NomActPosFile,'%d');

if ~silent 
    fprintf('...Alignment angle: %d\n',align_ang)
    fprintf('...Elevation angle: %d\n',RES.ElevAng)
end

%% Print nominal transforms

% Update M1 to FOS to include derotation by ElevAng
T_M1_FOS = vlCsMult(GRAVRES.T_M1A_FOSA,vlCsRotZ(-RES.ElevAng));
T_M1_FOS = GRAVRES.T_M1A_FOSA;

if ~silent
    fprintf('\nTransforms in Euler XYZ angles (displacement X,Y,Z followed by rotation Rx, then Ry, then Rz\n');
    fprintf('Linear dimensions are meters, angles are degrees\n');

    fprintf('\nNominal Transforms:\n')
    fprintf('Nominal M1 to M2: %.3f %.3f %.3f %.1f %.1f %.1f\n', vlCsEulerXYZ(GRAVRES.T_M1A_M2A));
    fprintf('Nominal M1 to M3: %.3f %.3f %.3f %.1f %.1f %.1f\n', vlCsEulerXYZ(GRAVRES.T_M1A_M3A));
    fprintf('Nominal M1 Aligned to FOS Aligned: %.3f %.3f %.3f %.1f %.1f %.1f\n', vlCsEulerXYZ(GRAVRES.T_M1A_FOSA))    
%     fprintf('Nominal M1 to FOS: %.3f %.3f %.3f %.1f %.1f %.1f\n', vlCsEulerXYZ(T_M1_FOS))
end

%% test of FOS transform 
% at 30 degrees align_ang this should result in Y = -1 once rotated into
% instrument coordinate system.
% and it does! (uncomment to check)

%fprintf('Test Transform:')
% GRAVRES.T_FOSDB_FOSDP = vlCsInvEulerXYZ([-1/2 -3^.5/2 0 0 0 RES.ElevAng-align_ang]);
% GRAVRES.T_FOSDB_FOSDP


%% Definitions
% Basis - at the aligned angle (usually 30 degrees)
% Perturbed - at the pointing zenith angle (0 to 65 or 90 degrees)

% Original - perfect alignment for perfectly placed primary mirror
% Aligned - perfect alignment to primary mirror best fit position
% Displaced - displacement under a load case before alignment takes place

% Notes on Alignment
% Focal surface B and P have with Y pointing towards zenith

%% Focal surface transformation derotated to zenith angles

% GRAVRES.T_FOSDB_FOSDP is the transform from displaced at alignment angle
% (30 deg for example) to displaced at perturbed angle (65 deg for example)
% Y axis is pointing along the M1 optical axis for both cases.
if ~silent
    fprintf('GRAVRES.T_FOSDB_FOSDP: %.3f %.3f %.3f %.1f %.1f %.1f\n', vlCsEulerXYZ(GRAVRES.T_FOSDB_FOSDP));
end

% remove the rotation from the transformation so that we are left with the
% displacement in the aligned coordinate system
FOS_EXYZ_B_P = vlCsEulerXYZ(GRAVRES.T_FOSDB_FOSDP) - [0 0 0 0 0 RES.ElevAng-align_ang];

if ~silent
    fprintf('FOS_EXYZ_B_P: %.3f %.3f %.3f %.1f %.1f %.1f\n', FOS_EXYZ_B_P);
end

% FOS_EXYZ_B_P
% Now rotate the displacement so that the FOS CS has Y pointing towards
% zenith, and take the inverse so that it is the transform of P to B

% to derotate a coordinate system transform:
% vlCsEulerXYZ(vlCsMult(vlCsRotZ(45),vlCsTrans([0,-1,0]'),vlCsRotZ(-45)))
% the first angle (+45) is the opposite sign of the desired new CS
% orientation.

FOS_Zenith_B_P = vlCsMult(vlCsRotZ(align_ang),vlCsInvEulerXYZ(FOS_EXYZ_B_P),vlCsRotZ(-align_ang));
FOS_Zenith_P_B = vlCsInv(FOS_Zenith_B_P);

if ~silent
    fprintf('FOS_Zenith_B_P: %.3f %.3f %.3f %.1f %.1f %.1f\n', vlCsEulerXYZ(FOS_Zenith_B_P));
end

%% Transform from Aligned M1 to the displaced M3
% Shoudln't, but need to remove rotation (need to check MFR code)
EXYZ_M1AP_M3DP = vlCsEulerXYZ(GRAVRES.T_M1AP_M3DP) - [0 0 0 RES.ElevAng-align_ang 0 0];
T_M1AP_M3DP = vlCsInvEulerXYZ(EXYZ_M1AP_M3DP);

%% Transform from Aligned M1 to the displaced M2
T_M1AP_M2DP = GRAVRES.T_M1AP_M2DP;

%% Transform from Aligned M1 to the displaced FOS CS
% Want the transform from M1 at the current telescope elevation angle to
% the displaced FOS.
% For example, basis elevation angle is 30 degrees, perturbed is 65 degrees

if ~silent
    fprintf('Transform from Aligned M1 to the displaced FOS CS\n');
end
% the following transform assumes that M1 is at the origin, whereas it
% really rotates about the elevation axis, 3.5 m way in Z.  Therefore we
% need to redo this transform with the proper rotations.

% There's a bug in the elevation angle difference in GRAVRES.T_M1AB_M1AP.
% Fixed here by zeroing out and then re-doing in the transform.
E_M1_B_P = vlCsEulerXYZ(GRAVRES.T_M1AB_M1AP);
E_M1_B_P(4) = 0;

if ~silent
    fprintf('GRAVRES.T_M1AB_M1AP, %.3f %.3f %.3f %.1f %.1f %.1f\n', vlCsEulerXYZ(GRAVRES.T_M1AB_M1AP))
    fprintf('E_M1_B_P: %.3f %.3f %.3f %.1f %.1f %.1f\n', E_M1_B_P);
    fprintf('ElevAng, AlignAng = %f,%f\n',RES.ElevAng, align_ang);
end

T_M1_B_P = vlCsInvEulerXYZ(E_M1_B_P);

T_M1AB_M1AP = vlCsMult(T_M1_B_P, vlCsTrans([0,0,3.5]'),vlCsRotX(RES.ElevAng-align_ang),vlCsTrans([0,0,-3.5]'));

if ~silent
    fprintf('T_M1AB_M1AP: %.3f %.3f %.3f %.1f %.1f %.1f\n', vlCsEulerXYZ(T_M1AB_M1AP));
end

% following line from the original code is in error (no need to add
% rotation)
% E_M1B_M1P = vlCsEulerXYZ(GRAVRES.T_M1AB_M1AP) + [0 0 0 RES.ElevAng-align_ang 0 0];

E_M1AB_M1AP = vlCsEulerXYZ(T_M1AB_M1AP);
T_M1AP_M1AB = vlCsInv(T_M1AB_M1AP);

if ~silent
    fprintf('T_M1AP_M1AB: %.3f %.3f %.3f %.1f %.1f %.1f\n', vlCsEulerXYZ(T_M1AP_M1AB));
end
%fprintf('T_M1AP_M1AB: %.3f %.3f %.3f %.1f %.1f %.1f\n', vlCsEulerXYZ(T_M1AP_M1AB));

% create the transform from M1 perfect at align_ang to M1 perfect at zenith
% remember that the rotation is about the elevation axis, so translate to
% the elevation axis, rotate, then translate back
T_M1BO_M1ZenithO = vlCsMult(vlCsTrans([0,0,3.5]'),vlCsRotX(-align_ang),vlCsTrans([0,0,-3.5]'));

T_M1BO_FOSBO = vlCsMult(T_M1BO_M1ZenithO, T_M1_FOS);

if ~silent
    fprintf('T_M1BO_FOSBO: %.3f %.3f %.3f %.1f %.1f %.1f\n', vlCsEulerXYZ(T_M1BO_FOSBO));
end

% Transform from M1P to FOSP
T_M1AP_FOSDP = vlCsMult(T_M1AP_M1AB, T_M1BO_FOSBO, FOS_Zenith_B_P);

if ~silent
    fprintf('T_M1AP_FOSDP: %.3f %.3f %.3f %.1f %.1f %.1f\n', vlCsEulerXYZ(T_M1AP_FOSDP));
end

if ~silent
    fprintf('Zenith Ang = %f\n', RES.ElevAng)
    fprintf('align Ang = %f\n', align_ang)
end

% fprintf('M1AB_M1AP = %f,%f,%f,%f,%f,%f\n',GRAVRES.EXYZ_M1AB_M1AP)
% fprintf('T_M1B_M1P\n');
% vlCsEulerXYZ(T_M1B_M1P)
% fprintf('T_M1_FOS\n');
% vlCsEulerXYZ(T_M1_FOS)
% fprintf('FOS_Zenith_P_B\n');
% vlCsEulerXYZ(FOS_Zenith_P_B)

% T_M1B_M1P
% T_M1_FOS
% FOS_Zenith_P_B

%% Results %%

% Print displaced transforms
if ~silent
    fprintf('\nResults:\n')
    fprintf('Transform from M1 Aligned to M2 Displaced: ')
    fprintf('%.4f %.4f %.4f %.2f %.2f %.2f\n',vlCsEulerXYZ(T_M1AP_M2DP))

    fprintf('Transform from M1 Aligned to M3 Displaced:')
    fprintf('%.4f %.4f %.4f %.2f %.2f %.2f\n',vlCsEulerXYZ(T_M1AP_M3DP))

    fprintf('Transform from M1 Aligned Perturbed to Focal Surface Displaced Perturbed: ')
    fprintf('%.4f %.4f %.4f %.2f %.2f %.2f\n', vlCsEulerXYZ(T_M1AP_FOSDP))
end

%% Return Values
T_M1M2 = T_M1AP_M2DP;
T_M1M3 = T_M1AP_M3DP;
T_M1FOS = T_M1AP_FOSDP;

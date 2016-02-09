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

if ~silent
    fprintf('\nTransforms in Euler XYZ angles (displacement X,Y,Z followed by rotation Rx, then Ry, then Rz\n');
    fprintf('Linear dimensions are meters, angles are degrees\n');

    fprintf('\nNominal Transforms:\n')
    fprintf('Nominal M1 to M2: %.3f %.3f %.3f %.1f %.1f %.1f\n', vlCsEulerXYZ(GRAVRES.T_M1A_M2A));
    fprintf('Nominal M1 to M3: %.3f %.3f %.3f %.1f %.1f %.1f\n', vlCsEulerXYZ(GRAVRES.T_M1A_M3A));

    fprintf('Nominal M1 to FOS: %.3f %.3f %.3f %.1f %.1f %.1f\n', vlCsEulerXYZ(T_M1_FOS))
end

%% test of FOS transform 
% at 30 degrees align_ang this should result in Y = -1 once rotated into
% instrument coordinate system.
% and it does! (uncomment to check)

%fprintf('Test Transform:')
% GRAVRES.T_FOSDB_FOSDP = vlCsInvEulerXYZ([-1/2 -3^.5/2 0 0 0 RES.ElevAng-align_ang]);
% GRAVRES.T_FOSDB_FOSDP

%% Focal surface transformation derotated to zenith angles
FOS_EXYZ_B_P = vlCsEulerXYZ(GRAVRES.T_FOSDB_FOSDP) - [0 0 0 0 0 RES.ElevAng-align_ang];
% FOS_EXYZ_B_P
FOS_Zenith_B_P = vlCsMult(vlCsRotZ(align_ang),vlCsInvEulerXYZ(FOS_EXYZ_B_P),vlCsRotZ(-align_ang));
FOS_Zenith_P_B = vlCsInv(FOS_Zenith_B_P);
FOS_Zenith_EXYZ_B_P = vlCsEulerXYZ(FOS_Zenith_B_P);

% to derotate a coordinate system transform:
% vlCsEulerXYZ(vlCsMult(vlCsRotZ(45),vlCsTrans([0,-1,0]'),vlCsRotZ(-45)))
% the first angle (+45) is the opposite sign of the desired new CS
% orientation.

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

E_M1B_M1P = vlCsEulerXYZ(GRAVRES.T_M1AB_M1AP) - [0 0 0 RES.ElevAng-align_ang 0 0];
T_M1B_M1P = vlCsInvEulerXYZ(E_M1B_M1P);
% vlCsEulerXYZ(T_M1B_M1P)

% Transform from M1P to FOSP
T_M1P_FOSP = vlCsMult(FOS_Zenith_P_B, T_M1_FOS, T_M1B_M1P);

%% Results %%

% Print displaced transforms
if ~silent
    fprintf('\nResults:\n')
    fprintf('Transform from M1 Aligned to M2 Displaced: ')
    fprintf('%.4f %.4f %.4f %.2f %.2f %.2f\n',vlCsEulerXYZ(T_M1AP_M2DP))

    fprintf('Transform from M1 Aligned to M3 Displaced:')
    fprintf('%.4f %.4f %.4f %.2f %.2f %.2f\n',vlCsEulerXYZ(T_M1AP_M3DP))

    fprintf('Transform from M1 Aligned to Focal Surface Displaced: ')
    fprintf('%.4f %.4f %.4f %.2f %.2f %.2f\n', vlCsEulerXYZ(T_M1P_FOSP))
end

%% Return Values
T_M1M2 = T_M1AP_M2DP;
T_M1M3 = T_M1AP_M3DP;
T_M1FOS = T_M1P_FOSP;

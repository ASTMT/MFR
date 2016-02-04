% File: nrSciCalcResults.m
%
% Syntax: nrSciCalcResults
%
% Description:
%       Calculates results in RES data structure for the Merit Function
%       (nrSciMerit). 
%
% Input Parameters:
%       None.
%
% Output Parameters:
%       None
%
% Required Global Data Structures:
%       OC
%       IM
%       MA
%       NRCIM
%       CST - The CST structure being validated must be loaded.
%       OC  - The OC data structure is required.
%       RES - The Merit Function result data structure.  This is
%       initialized in vlSciMerit
%
% Required Data Files:
%       N/A
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%       This function is a rewritten version of vlSciCalcResults using the
%       NRCIM data structures

%
% Revision History
%
% $Id: 
%
% INDENT-OFF*
%
% INDENT-ON*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%           Herzberg Institute of Astrophysics                  %%%%%
%%%%%%      Astronomy Technology Research Group - Victoria           %%%%%
%
% (c) <2009>				        (c) <2009>
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

function nrSciCalcResults

global RES
global OC
global IM
global CST
global MA
global NRCIM
global TMTCS

%% Calculate science merit function results

fprintf('....Calculating Science Merit Function Results\n');

%% Set up some constants and definitions

asecPerRad = 180*3600/(pi);
radPerDeg = pi/180;

% identify the surfaces rather than using integer numbers
M1_surfnum = 1;
M2_surfnum = 2;
M3_surfnum = 3;
foc_surfnum = 4;

% Now find the CST index row number for the sequential surfaces
M2surfidx = find(CST.SurfList==M2_surfnum);
M3surfidx = find(CST.SurfList==M3_surfnum);
focsurfidx = find(CST.SurfList==foc_surfnum);

NumSegs = NRCIM.SegsPerSurf(M1_surfnum);
M1cstlist = CST.SurfBool(:,M1_surfnum);

%%
CSTZPZ = vlCsInv(RES.CSTZZP);

ActStroke = RES.ActStroke;
RES.ACT_MAX = max(abs(ActStroke));
% Write actuator values to a file

% save(RES.ActFileName,'ActStroke','-ASCII');
RES.ACT_RMS = norm(ActStroke)/sqrt(length(ActStroke));
RES.ACT_MEAN = mean(ActStroke);

%%  Calculate the uncorrected image quality due to the primary mirror
%  decenter and rotation only IQ_PRIM_*.  I.E. take the least squares fit
%  primary and calculate the image quality with no error on any of the
%  sequential mirrors and with tip/tilt piston removed from the primary
%  segments.


%% calculate the sp array
% ActPosAlignY only includes primary mirror motions, so there is no need to
% zero out the perturbations for other surfaces

Y = RES.ActPosAlignY;
[sp TS] = nrCimCalcPerturbations(Y);

% If the segment SP basis position vector was supplied using the
% '-NomAct' flag, then subtract this from the sp vector.
if ~isempty(RES.NomActPosFile)
    sp = sp - RES.ActPosNom.SP;
else % save the SP vector used above
    RES.ActPosNom.SP = sp;
end

RES.SPBeforeAlignment = sp;

% Define the following results giving the TMT Segment Coordinate system
% perturbations for input to the JPL M1CS analysis 
    % RES.CST_M1CS_TMTM1Seg_Original
    % RES.CST_M1CS_TMTM1Seg_Perturbed
    % RES.EulerXYZ_TMTM1Seg_orig_to_pert
    
nrSciGetTMTM1CSPerturbations;

% Create sp vectors for the following cases:
%   - only primary mirror segment decenter and rotation (decrot)
%       In this case Dz = 0 and Rx and Ry are set for the Dx,Dy position
%   - only primary mirror segment decenter (dec)
%       In this case Dz = Rz = 0 and Rx and Ry are set for the Dx,Dy position
%   - only primary mirror segment rotation (rot)
%       In this case Dx = Dy = Dz = Rx = Ry = 0 and Rz is non-zero

% sp column order is Dx,Dy,Dz,Rx,Ry,Rz
sp_decrot = sp;
% zero out the Dz and Rx,Ry
% Dx and Dy are inputs below and Rz is preserved
%sp_decrot(:,3:5) = 0;

% compensate segment decenter with tip-tilt

for ii = 1:CST.Num
    % Does the CST belong to the primary mirror?
    if CST.Meta(ii).NrcimSurf == M1_surfnum
        % set the tip-tilt for the new position
        %[sp_decrot(ii,4), sp_decrot(ii,5), Z, Set] =
        %vlOpDecComp(CST.ZM(:,:,ii),sp_decrot(ii,:));
        sp_decrot(ii,:) = vlOpDecComp(CST.ZM(:,:,ii),sp_decrot(ii,:));
    end
end

% save the compensated primary mirror sp values.  This includes rows for
% the sequential surfaces, but they are zeroed out
RES.sp_compensated_M1 = sp_decrot(logical(CST.SurfBool(:,M1_surfnum)),:);

% sp_decrot
sp_dec = sp_decrot;
sp_dec(:,6) = 0;

sp_rot = sp_decrot;
sp_rot(:,1:5) = 0;

% vectorize the sp nx6 arrays
vec_sp_decrot = sp_decrot';
vec_sp_decrot = vec_sp_decrot(:);

vec_sp_dec = sp_dec';
vec_sp_dec = vec_sp_dec(:);

vec_sp_rot = sp_rot';
vec_sp_rot = vec_sp_rot(:);

% calculate the OPD 
% want only the contribution of the telescope, so subtract the OPDzero

RES.opd_decrot = vlLmRunLom(vec_sp_decrot) - OC.LOM.OPDzero;
RES.opd_dec = vlLmRunLom(vec_sp_dec) - OC.LOM.OPDzero; 
RES.opd_rot = vlLmRunLom(vec_sp_rot) - OC.LOM.OPDzero;

% calculate the RMS wavefront
RES.IQ_DECROT_RMS = vlOpWtoRMS(RES.opd_decrot);
RES.IQ_DEC_RMS = vlOpWtoRMS(RES.opd_dec);
RES.IQ_ROT_RMS = vlOpWtoRMS(RES.opd_rot);

if RES.DebugMode
    % Calculate expected nm RMS from Mast's equations
    Mast_decrms = 0.707*0.595*(OC.SegSize/2)^2*abs(OC.PrimCC)*OC.EPD/2*RES.PRIM_DEC_RMS/OC.PrimRadCurv^3;
    Mast_rotrms = 1/sqrt(3)*0.341*abs(OC.PrimCC)*(OC.SegSize/2)^2*(OC.EPD/2)^2*RES.PRIM_ROT_RMS*radPerDeg/OC.PrimRadCurv^3;
    fprintf('Check from Terry Masts equations: RMS decenter = %.3e, RMS rotation = %.3e\n',Mast_decrms, Mast_rotrms);
end

%% Writes SP vector to file & Plots OPD

if RES.OutputSP
   save([RES.outFileRoot '_sp_decrot.txt'],'sp_decrot','-ASCII');
   save([RES.outFileRoot '_sp_dec.txt'],'sp_dec','-ASCII');
   save([RES.outFileRoot '_sp_rot.txt'],'sp_rot','-ASCII');
end

RES.PRIM_ROT_MAX = max(abs(RES.sp_compensated_M1(:,6)));
RES.PRIM_ROT_MEAN = mean(RES.sp_compensated_M1(:,6));
RES.PRIM_ROT_RMS = norm(RES.sp_compensated_M1(:,6))/sqrt(NumSegs);
dec = sqrt(RES.sp_compensated_M1(:,1).^2+RES.sp_compensated_M1(:,2).^2);
RES.PRIM_DEC_MAX = max(abs(dec));
RES.PRIM_DEC_MEAN = mean(dec);
RES.PRIM_DEC_RMS = norm(dec)/sqrt(NRCIM.Surf(1).ZObjNum);
RES.PRIM_DEC_X_MAX = max(abs(RES.sp_compensated_M1(:,1)));
RES.PRIM_DEC_Y_MAX = max(abs(RES.sp_compensated_M1(:,2)));
RES.PRIM_DEC_X_MEAN = mean(RES.sp_compensated_M1(:,1));
RES.PRIM_DEC_Y_MEAN = mean(RES.sp_compensated_M1(:,2));
RES.PRIM_DEC_X_RMS = norm(RES.sp_compensated_M1(:,1))/sqrt(NumSegs);
RES.PRIM_DEC_Y_RMS = norm(RES.sp_compensated_M1(:,2))/sqrt(NumSegs);


%% Generate the TMT Coordinate Systems
% % call vlCsTmtCST to define the TMT coordinate systems.  Let the routine
% % default to Az=0, El=0, and 'Northern' hemisphere since we don't know
% % these values. Note that these inputs only change coordinate transforms 
% % below the primary mirror which we don't need for the merit function.
% 
% if vlSciMeritGetVersion >= 6.0
%     distances = eval(IM.GetDistancesFile); 
% else
%     distances = vlCsM1pGetDistances_3mirr; 
% end
% TmtCs = vlCsTmtCst(distances);

% want to find the transforms between Z and TMT M1, M2, M3 and Focus. The
% final results should be reported in these coordinate systems.

RES.CST_Z_TMTM1 = TMTCS.ZM1;
RES.CST_Z_TMTM2 = vlCsMult(TMTCS.ZM1,TMTCS.M1M2);
RES.CST_Z_TMTM3 = vlCsMult(TMTCS.ZM1,TMTCS.M1M3);
RES.CST_Z_TMTFOS = vlCsMult(TMTCS.ZM1,TMTCS.M1FOC);

%% if not a NomAct solution, calculate the results needed for the single case analysis
if ~RES.NomAct
    % Calculate how M1 has moved relative to the undisplaced position

    % Calculate the transform between Tmt M1 and Tmt M1'
    [RES.CSTM1P,RES.EXYZM1P] = vlCsRigidFrameTrans([eye(3) zeros(3,1)],RES.CSTZZP,TMTCS.ZM1);

    % Calculate the Center of Curvature of displaced primary in the undisplaced TMT M1 CS
    RES.COC = vlCsPMult(vlCsMult(vlCsInv(TMTCS.ZM1),RES.CSTZZP),[0 0 -OC.PrimRadCurv]');
    RES.COC_X = RES.COC(1);
    RES.COC_Y = RES.COC(2);
    RES.COC_Z = RES.COC(3);

    % Calculate how M2 has moved relative to the Z prime coordinate system, and
    % what correction is required to re-align it with the optical system based
    % on the new M1 position.

    % Coordinate frame definitions
    % M2O - NRCIM M2 CST Original Position
    % M2A - NRCIM M2 CST Alignment Position
    % M2D - NRCIM Displaced Position of M2
    % TMTM2O - TMT M2 CST Original Position
    % TMTM2A - TMT M2 CST Alignment Position
    % TMTM2D - TMT M2 CST Displaced Position

    % Calculate the transform between the NRCIM M2 and TMT M2 frames
    CST_M2_TMTM2 = vlCsMult(vlCsInv(CST.ZM(:,:,M2surfidx)),TMTCS.ZM1,TMTCS.M1M2);

    % Calculate the transform between Z, Z Prime and the M2 Displaced Position
    CSTZM2D = vlCsMult(CST.ZC(:,:,M2surfidx),vlCsInv(RES.TSFULL.CPC(:,:,M2surfidx)),CST.CM(:,:,M2surfidx));
    CSTZPM2D = vlCsMult(CSTZPZ,CSTZM2D);

    % We want the transform from ZP to M2 Aligned to be the ideal position 
    CSTZPM2A = CST.ZM(:,:,M2surfidx);

    % now we have the ideal M2 CS and the displaced M2 CS relative to the
    % Z prime coordinate system. These M2 CS's are the integrated model 
    % definitions, but we want the answer in TMT CS's.

    % Want the transform from M2O to M2D
    CST_M2O_M2D = vlCsMult(vlCsInv(CSTZPM2A),CSTZM2D);
    %RES.EULM2N = vlCsEulerXYZ(CST_M2O_M2D);
    RES.EULM2O_M2D = vlCsEulerXYZ(CST_M2O_M2D);

    % Get the transform from M2A to M2D
    CST_M2A_M2D = vlCsMult(vlCsInv(CSTZPM2A),CSTZPM2D);
    %RES.EULM2P = vlCsEulerXYZ(CST_M2A_M2D);
    RES.EULM2A_M2D = vlCsEulerXYZ(CST_M2A_M2D);

    % Get the transform from M1O to M2D
    CST_M1O_M2D = vlCsMult(TMTCS.M1M2,CST_M2O_M2D);
    RES.EULM1O_M2D = vlCsEulerXYZ(CST_M1O_M2D);

    % Get the transform from M1A to M2D
    CST_M1A_M2D = vlCsMult(TMTCS.M1M2,CST_M2A_M2D);
    RES.EULM1A_M2D = vlCsEulerXYZ(CST_M1A_M2D);

    % Calculate how M3 has moved relative to the Z prime coordinate system, and
    % what correction is required to re-align it with the optical system based
    % on the new M1 and M3 positions.

    % M3O - M2 CST Original Position
    % M3A - M2 CST Alignment Position
    % M3D - Displaced Position of M2

    % Calculate the transform between Z, Z Prime and M3 Displaced Position
    CSTZM3D = vlCsMult(CST.ZC(:,:,M3surfidx),vlCsInv(RES.TSFULL.CPC(:,:,M3surfidx)),CST.CM(:,:,M3surfidx));
    CSTZPM3D = vlCsMult(CSTZPZ,CSTZM3D);

    % We want the transform from ZP to M3 Aligned to be the ideal position 
    CSTZPM3A = CST.ZM(:,:,M3surfidx);

    % now we have the ideal M3 CS and the displaced M3 CS relative to the
    % Z prime coordinate system. These M3 CS's are the integrated model 
    % definitions, but we want the answer in TMT CS's.
    % The outputs of the following call give the gross motion of M3 relative to
    % the displaced M1 coordinate system (Z Prime or M1').

    % Want the transform from M3O to M3D
    CST_M3O_M3D = vlCsMult(vlCsInv(CSTZPM3A),CSTZM3D);
    %RES.EULM2N = vlCsEulerXYZ(CST_M2O_M2D);
    RES.EULM3O_M3D = vlCsEulerXYZ(CST_M3O_M3D);

    % Get the transform from M3A to M3D
    CST_M3A_M3D = vlCsMult(vlCsInv(CSTZPM3A),CSTZPM3D);
    %RES.EULM2P = vlCsEulerXYZ(CST_M2A_M2D);
    RES.EULM3A_M3D = vlCsEulerXYZ(CST_M3A_M3D);

    % Get the transform from M1O to M3D
    CST_M1O_M3D = vlCsMult(TMTCS.M1M3,CST_M3O_M3D);
    RES.EULM1O_M3D = vlCsEulerXYZ(CST_M1O_M3D);

    % Get the transform from M1A to M3D
    CST_M1A_M3D = vlCsMult(TMTCS.M1M3,CST_M3A_M3D);
    RES.EULM1A_M3D = vlCsEulerXYZ(CST_M1A_M3D);

    % Calculate how the focal plane has moved relative to the Z prime coordinate system, and
    % what correction is required to re-align it with the optical system based
    % on the new M1 and M2 positions.

    % Calculate the transfrom between Z, Z prime and the focal plane Displaced
    % Positions
    CSTZFOSD = vlCsMult(CST.ZC(:,:,focsurfidx),vlCsInv(RES.TSFULL.CPC(:,:,focsurfidx)),CST.CM(:,:,focsurfidx));
    CSTZPFOSD = vlCsMult(CSTZPZ,CSTZFOSD);

    % We want the transform from ZP to the focal plane Aligned to be the ideal position 
    CSTZPFOSA = CST.ZM(:,:,focsurfidx);

    % now we have the ideal FOS CS and the displaced FOS CS relative to the
    % Z prime coordinate system. These FOS CS's are the integrated model 
    % definitions, but we want the answer in TMT CS's.
    % The outputs of the following call give the gross motion of FOS relative to
    % the displaced M1 coordinate system (Z Prime or M1').
    % Want the transform from FOSO to FOSD
    CST_FOSO_FOSD = vlCsMult(vlCsInv(CSTZPFOSA),CSTZFOSD);
    %RES.EULM2N = vlCsEulerXYZ(CST_M2O_M2D);
    RES.EULFOSO_FOSD = vlCsEulerXYZ(CST_FOSO_FOSD);

    % Get the transform from FOSA to FOSD
    CST_FOSA_FOSD = vlCsMult(vlCsInv(CSTZPFOSA),CSTZPFOSD);
    %RES.EULM2P = vlCsEulerXYZ(CST_M2A_M2D);
    RES.EULFOSA_FOSD = vlCsEulerXYZ(CST_FOSA_FOSD);

    % Get the transform from M1O to FOSD
    CST_M1O_FOSD = vlCsMult(TMTCS.M1FOC,CST_FOSO_FOSD);
    RES.EULM1O_FOSD = vlCsEulerXYZ(CST_M1O_FOSD);

    % Get the transform from M1A to FOSD
    CST_M1A_FOSD = vlCsMult(TMTCS.M1FOC,CST_FOSA_FOSD);
    RES.EULM1A_FOSD = vlCsEulerXYZ(CST_M1A_FOSD);

    % [CST_ZP_TMTFOS_TMTFOSP,EXYZ_ZP_TMTFOS_TMTFOSP] = vlCsRigidFrameTrans(CSTZPFOS,CSTZPFOSP,RES.CST_Z_TMTFOS);
    % 
    % % now calculate the gross motion from the original Z coordinate system
    % [CST_Z_TMTFOS_TMTFOSP,EXYZ_Z_TMTFOS_TMTFOSP] = vlCsRigidFrameTrans(CSTZPFOS,CSTZFOSP,RES.CST_Z_TMTFOS);
    % 
    % % RES.M3P gives euler angles of the position of M3 relative to the aligned primary mirror
    % RES.EULFOSP = EXYZ_ZP_TMTFOS_TMTFOSP;
    % 
    % % RES.M2P gives euler angles of the position of M2 relative to the original primary mirror
    % % position
    % RES.EULFOSN = EXYZ_Z_TMTFOS_TMTFOSP;
end
%% Get the node numbers and positions that correspond to the actuator
% position vector.  The actuator vector is in the order of the CSTs for the
% primary mirror surface

% get the list of CST positions for M1, in order
RES.ActNodeNumList = zeros(NumSegs*3,1);    % Ansys node numbers

% generate n x 3 list

RES.ActNodeNumList = [CST.Meta(logical(M1cstlist)).NodeNums];
% change to n*3 x 1 list
RES.ActNodeNumList = RES.ActNodeNumList(:);


% ActNodePosList is a NumSegs*3 x 3 array of actuator positions in the Z
% coordinate system.
% ActNodePosListTMTM1 is ActNodePosList in TMT M1 coordinates

RES.ActNodePosList = zeros(length(RES.ActNodeNumList),3);
RES.ActNodePosListTMTM1 = zeros(length(RES.ActNodeNumList),3);

for ii = 1:length(RES.ActNodeNumList)
    RES.ActNodePosList(ii,1:3) = vlCsPMult(CST.ZG,MA.NodePosG(RES.ActNodeNumList(ii),:)');
    P = vlCsPMult(vlCsInv(RES.CST_Z_TMTM1),RES.ActNodePosList(ii,:)');
    RES.ActNodePosListTMTM1(ii,:) = P';
end
                                                                
%% Align telescope and correct LOS image and focus          

% Alignment plan
% - Telescope distorts, primary and secondary are re-aligned
% - Tertiary and Focal plane sag are added
% - Focal plane is transformed based on primary mirror alignment
% - Tertiary is then tilted to repoint at the center of the focal plane


% define the rows in the SP vector for sequential surfaces
SEC_START = (M2surfidx-1)*6 + 1;
SEC_PISTON = SEC_START + 2;

TER_START = (M3surfidx-1)*6 + 1;
TER_X_TILT_POS = TER_START + 3;
TER_Y_TILT_POS = TER_START + 4;

FOS_START = (focsurfidx-1)*6 + 1;

%% Calculate ZER20 sensitivity to secondary piston

piston = 0.0001; 

% Calculate the focus term before and after piston
% Generate the sp vector
sp_OSM = zeros(CST.Num*6,1);
opd = vlLmRunLom(sp_OSM);
zernikes = vlOpWtoZer(opd);
focus = zernikes(4);

sp_OSM(SEC_PISTON) = piston;
% calculate the opd and zernikes
opd = vlLmRunLom(sp_OSM);
zernikes = vlOpWtoZer(opd);
focus = zernikes(4)-focus;

% calculate the sensitivity
RES.SEN_SEC_TZ_ZER20 = focus/piston;  % zer per m

%% Calculate sensitivities to tertiary mirror tilt
% NOTE: THE FOLLOWING CODE IS COMMENTED OUT UNTIL IT CAN BE DEBUGGED AND
% VALIDATED

% tilt = 1/60; % 1 arcminute
% 
% %
% % first tilt about tertiary X axis
% %
% % Generate the sp vector
% sp_OSM = zeros(CST.Num*6,1);
% sp_OSM(TER_X_TILT_POS) = tilt;
% % calculate the opd and zernikes
% opd = vlLmRunLom(sp_OSM);
% zernikes = vlOpWtoZer(opd);
% x = zernikes(2)*OC.FL/OC.EPD*4;
% y = zernikes(3)*OC.FL/OC.EPD*4;
% % Calculate the sensitivities of TER rotation to FP translation
% RES.SEN_TER_RX_FX = x/tilt;   % result is in m/deg
% RES.SEN_TER_RX_FY = y/tilt;
% 
% %
% % then tilt about tertiary Y axis
% %
% sp_OSM = zeros(CST.Num*6,1);
% sp_OSM(TER_Y_TILT_POS) = tilt;
% % calculate the opd and zernikes
% opd = vlLmRunLom(sp_OSM);
% zernikes = vlOpWtoZer(opd);
% x = zernikes(2)*OC.FL/OC.EPD*4;
% y = zernikes(3)*OC.FL/OC.EPD*4;
% % Calculate the sensitivities of TER rotation to FP translation
% RES.SEN_TER_RY_FX = x/tilt;   % result is in m/deg
% RES.SEN_TER_RY_FY = y/tilt;
% 
% %% Create the SP vector for the telescope without tertiary correction
% 
% % start with the SP vector of primary segment perturbations
% sp = vec_sp_decrot;
% 
% % add the secondary motions relative to the alignment coordinates. 
% sp(SEC_START:SEC_START+5) = RES.EULM2A_M2D;
% 
% % add the tertiary motions relative to the alignment coordinates. 
% sp(TER_START:TER_START+5) = RES.EULM3A_M3D;
% 
% % add the focal plane motions relative to the alignment coordinates. 
% sp(FOS_START:FOS_START+5) = RES.EULFOSA_FOSD;
% 
% %% find LOS centroid of image
% 
% 
% % Now calculate the image position
% % calculate the opd and zernikes
% 
% opd = vlLmRunLom(sp);
% zernikes = vlOpWtoZer(opd);
% RES.LOS_X = -zernikes(2)*OC.FL/OC.EPD*4;
% RES.LOS_Y = -zernikes(3)*OC.FL/OC.EPD*4;
% 
% %% find LOS centroid of image
% 
% 
% % Write 2 x 2 matrix that converts required image motions into tertiary
% % rotations
% FocPos_to_TerRot = [RES.SEN_TER_RX_FX RES.SEN_TER_RY_FX; RES.SEN_TER_RX_FY RES.SEN_TER_RY_FY]^-1;
% TerRot = FocPos_to_TerRot*[RES.LOS_X;RES.LOS_Y];
% RES.TER_POINT_RX = TerRot(1);
% RES.TER_POINT_RY = TerRot(2);
% 
% 
% %% Find secondary piston to correct focus
% RES.SEC_FOCUS = -zernikes(4)/RES.SEN_SEC_TZ_ZER20;
% 
% %%
% sp(TER_X_TILT_POS) = sp(TER_X_TILT_POS)+RES.TER_POINT_RX;
% sp(TER_Y_TILT_POS) = sp(TER_Y_TILT_POS)+RES.TER_POINT_RY;
% sp(SEC_PISTON) = sp(SEC_PISTON) + RES.SEC_FOCUS;
% 
% %% Calculate the plate scale
% 
% %Get primary to secondary distance, D
% Dideal = abs(CST.ZM(3,4,M2surfidx));
% Dalign = Dideal + RES.SEC_FOCUS;
% 
% RES.PlateIdeal = (OC.PrimRadCurv/2)/(2*Dideal/OC.SecRadCurv - OC.PrimRadCurv/OC.SecRadCurv - 1)/asecPerRad;
% RES.PlateAlign = (OC.PrimRadCurv/2)/(2*Dalign/OC.SecRadCurv - OC.PrimRadCurv/OC.SecRadCurv - 1)/asecPerRad;
% 
% %% Now recalculate opd and save zernikes for the corrected case
% 
% RES.CorrOPD = vlLmRunLom(sp);
% 
% % calculate the encircled energy due to decenter and rotation only.
% [psf,icr] = vlOpWtoPSFFFT(RES.CorrOPD,OC.LOM.IAD,OC.Wave(OC.LOM.WaveNum),OC.FFT.Zp,OC.EPD);  
% 
%  % want ee in arcseconds diameter, so must divide by imagescale
% ee = vlOpPSFtoEE(psf,[0.5 0.8],icr)/OC.ImageScale;
% RES.IQ_ALL_50EE = ee(1); 
% RES.IQ_ALL_80EE = ee(2);
% 
% % calculate the Strehl Ratio
% RES.IQ_ALL_STREHL = vlOpPSFtoStrehl(psf);
% 
% % calculate the Zernikes
% zernikes = vlOpWtoZer(RES.CorrOPD);
% 
% % Calculate the x,y position for debugging
% x = zernikes(2)*OC.FL/OC.EPD*4;
% y = zernikes(3)*OC.FL/OC.EPD*4;
% 
% RES.ZERNIKES = zernikes;
% RES.ZER11C = zernikes(2);
% RES.ZER11S = zernikes(3);
% RES.ZER20 = zernikes(4);
% RES.ZER22S = zernikes(5);
% RES.ZER22C = zernikes(6);
% RES.ZER31S = zernikes(7);
% RES.ZER31C = zernikes(8);
% RES.ZER33S = zernikes(9);
% RES.ZER33C = zernikes(10);
% RES.ZER40 = zernikes(11);
% RES.ZER42S = zernikes(12);
% RES.ZER42C = zernikes(13);
% RES.ZER44C = zernikes(14);
% RES.ZER44S = zernikes(15);
% 
% RES.IQ_ALL_RMS = vlOpWtoRMS(RES.CorrOPD);
% 
% %% Now recalculate opd and save zernikes for the Un-corrected case
% 
% RES.UnCorrOPD = vlLmRunLom(RES.UncorrectedSp);
% 
% % calculate the encircled energy due to decenter and rotation only.
% [psf,icr] = vlOpWtoPSFFFT(RES.UnCorrOPD,OC.LOM.IAD,OC.Wave(OC.LOM.WaveNum),OC.FFT.Zp,OC.EPD);  
% 
%  % want ee in arcseconds diameter, so must divide by imagescale
% ee = vlOpPSFtoEE(psf,[0.5 0.8],icr)/OC.ImageScale;
% RES.IQ_UNCORR_50EE = ee(1); 
% RES.IQ_UNCORR_80EE = ee(2);
% 
% % calculate the Strehl Ratio
% RES.IQ_UNCORR_STREHL = vlOpPSFtoStrehl(psf);
% 
% % calculate the Zernikes
% zernikes = vlOpWtoZer(RES.UnCorrOPD);
% 
% % Calculate the x,y position for debugging
% x = zernikes(2)*OC.FL/OC.EPD*4;
% y = zernikes(3)*OC.FL/OC.EPD*4;
% 
% RES.UNCORR_ZERNIKES = zernikes;
% 
% RES.IQ_UNCORR_RMS = vlOpWtoRMS(RES.UnCorrOPD);

% We're done!
%fprintf('....done in %.1f Seconds\n\n',toc);

% end of nrSciCalcResults

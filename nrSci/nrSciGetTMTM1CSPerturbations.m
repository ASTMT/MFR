% File: <nrSciGetTMTM1CSPerturbations.m>
%
% Syntax: nrSciGetTMTM1CSPerturbations
%
% Description:
%       Calculates the M1 segment perturbations of the official TMT segment
%       coordinate systems written in M1CS coordinates.
%
% Input Parameters:
%       None
%
% Output Parameters:
%       None
%
% Required Global Data Structures:
%       RES
%       CST
%       TMTCS
%       NRCIM
%
% Required Data Files:
%       None
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: nrSciGetTMTM1CSPerturbations.m,v 1.2 2011/02/10 18:19:45 roberts Exp $
%
% INDENT-OFF*

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


function nrSciGetTMTM1CSPerturbations

global RES
global CST
global TMTCS
global NRCIM
% load(NRCIMmatname);
% load(RESFileName);

CST_M1CS_TMTSEG = TMTCS.M1Seg;

%% Definitions
M1surf = 1;         % Identifier for M1 segment CSTs in CST.SurfBool and CST.SurfList
NumSegs = NRCIM.Surf(M1surf).ZObjNum;

%% Calculate NRCIM segment coordinate systems in M1CS

% Calculate the uncorrected (no tip-tilt or piston correction)
% perturbations of the segments in M coordinates.
% Y = RES.ActPosAlignY;
% [sp TS] = nrCimCalcPerturbations(Y);

sp = RES.SPBeforeAlignment;

% Get all the M1 coordinate systems
CSTM1.ZM = CST.ZM(:,:,logical(CST.SurfBool(:,1)));
sp = sp(logical(CST.SurfBool(:,1)),:);

CSTM1.M1CSM = zeros(3,4,NumSegs);

for ii = 1:NumSegs
    CSTM1.M1CSM(:,:,ii) = vlCsMult(vlCsInv(TMTCS.ZM1),CSTM1.ZM(:,:,ii));
    CSTM1.M1CSMP(:,:,ii) = vlCsMult(vlCsInv(TMTCS.ZM1),CSTM1.ZM(:,:,ii),vlCsInvEulerXYZ(sp(ii,:)));
end

%% make a map of NRCIM to TMT segment coordinate systems

% The segment order in the MFR is the same as the official TMT segment
% coordinate system order.  Therefore this section of code is not required.

% matches = zeros(NumSegs,1);
% 
% for nrcimseg = 1:NumSegs
%     num_matches = 0;
%     for tmtseg = 1:NumSegs
%         if norm(CST_M1CS_TMTSEG(:,4,tmtseg)-CSTM1.M1CSM(:,4,nrcimseg)) < 0.1
%             num_matches = num_matches+1;
%             matches(nrcimseg) = tmtseg;
%         end
%     end
%     if num_matches ~= 1
%         error('num_matches not equal to 1');
%     end
% end

%% Find the TMT segment perturbed coordinate system 

% Calculate the M1CS to TMTSEG' transformation

% File: vlCsRigidFrameTrans.m, NRCIM Toolbox
%
% Syntax: [CSTBBP,EBBP] = vlCsRigidFrameTrans(CSTZA,CSTZAP,CSTZB)

TMTSEG_SP = zeros(NumSegs,6);

for ii = 1:NumSegs
    % CSTSSP is short for CST TMTSeg to TMTSeg'
    [CSTSSP,ESSP] = vlCsRigidFrameTrans(CSTM1.M1CSM(:,:,ii),CSTM1.M1CSMP(:,:,ii),CST_M1CS_TMTSEG(:,:,ii));
    CST_M1CS_TMTSEGP(:,:,ii) = vlCsMult(CST_M1CS_TMTSEG(:,:,ii),CSTSSP);
    TMTSEG_SP(ii,:) = vlCsEulerXYZ(vlCsMult(vlCsInv(CST_M1CS_TMTSEG(:,:,ii)),CST_M1CS_TMTSEGP(:,:,ii)));
end

RES.CST_M1CS_TMTM1Seg_Original = CST_M1CS_TMTSEG;
RES.CST_M1CS_TMTM1Seg_Perturbed = CST_M1CS_TMTSEGP;
RES.EulerXYZ_TMTM1Seg_orig_to_pert = TMTSEG_SP;

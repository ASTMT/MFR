% File: nrCimCalcPerturbations
%
% Syntax: [SP TS] = nrCimCalcPerturbations(Y)
%
% Description:
%       Post Processor for the telescope structural dynamics.
%       Finds sequential surface perturbations.
%
% Input Parameters:
%       Y         - (MA.r x 1) Structural Dynamics displacements in nodal
%                   coordinates.
%
% Output Parameters:
%       SP      -   [CST.Num,6] The cst perturbations equivalent to 
%                   the Euler angles to get from M to M' for each 
%                   non-sequential surface.  The row order in the SP vector
%                   is the same order as surfaces identified in the
%                   CST.Meta data structure.
%       TS      -   TS data structure containing the CST perturbations for
%                   this analysis.
%
% Required Global Data Structures:
%       CST
%       TS - The TS data structure is modified by this routine
%
% Required Data Files:
%       none
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
%
% Revision History
%
% $Id: nrCimCalcPerturbations.m,v 1.2 2010/10/30 18:43:56 roberts Exp $
%
% INDENT-OFF*
% $Log: nrCimCalcPerturbations.m,v $
% Revision 1.2  2010/10/30 18:43:56  roberts
% Fixed error in 1node code
%
% Revision 1.1  2008/04/30 19:35:45  roberts
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


function [SP TS] = nrCimCalcPerturbations(Y)

global CST
global TS

% The Y vector contains the perturbations in nodal coordinates

for CstIdx = 1:CST.Num
    switch CST.Meta(CstIdx).AType
        case 'Triad'
            % Get the points in the C coordinate system
            CP1 = nrCimCP(CstIdx,1);
            CP2 = nrCimCP(CstIdx,2);
            CP3 = nrCimCP(CstIdx,3);

            % Get the perturbed points in the C coordinate system
            CP1P = CP1 + nrCimND(CstIdx,1,Y);
            CP2P = CP2 + nrCimND(CstIdx,2,Y);
            CP3P = CP3 + nrCimND(CstIdx,3,Y);

            % Compute transform for points in C' to C.
            TransCCP = vlCsTriadTrans(CP1P,CP2P,CP3P);
            
        case '1Node'
            CP = nrCimND(CstIdx,1,Y);
            TransCCP = vlCsInvEulerXYZ(CP');
        otherwise 
            error('Bad AType');
    end
    
    TS.CPC(:,:,CstIdx) = vlCsInv(TransCCP);
    % Compute the transform matrix to give the following transforms:
    % M --> C --> C' --> M'
    TS.MPM(:,:,CstIdx) = vlCsMult(vlCsInv(CST.CM(:,:,CstIdx)), ...
                                TS.CPC(:,:,CstIdx), ...
                                CST.CM(:,:,CstIdx));

    % Compute the transform matrix to go from C' to Z
    TS.ZCP(:,:,CstIdx) = vlCsMult(CST.ZC(:,:,CstIdx), TransCCP);

    % Compute the Euler angles
    TS.SSP(CstIdx,:) = vlCsEulerXYZ(vlCsInv(TS.MPM(:,:,CstIdx)));
end

if size(TS.SSP) ~= [CST.Num 6]
    error('Invalid size for SSP');
end

SP = TS.SSP;

% end of nrCimCalcPerturbations

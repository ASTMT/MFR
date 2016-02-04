% File: <nrSciMeritPlotSegDisp.m>
%
% Syntax: [xyuv h] = nrSciMeritPlotSegDisp(pmsp, plotflag)
%
% Description:
%       Produces a vector plot of segment displacement on the primary
%       mirror in the M1CS coordinate system.
%
% Input Parameters:
%       pmsp      - (OC.NumNonSeqs x 8) The standard pmsp data structure
%                   including x and y decenter data.  
%       plotflag  - (binary) 1 if figure should be created, 0 if only xyuv
%                   is calculated and returned.
%
% Output Parameters:
%       xyuv      - (OC.NumNonSeqs x 8) A table of segment x,y position and
%                   u,v displacements in the Z coordinate system.  This is
%                   the data that is plotted in the figure.
%       h         - handle to the figure, set to 0 if plotflag = 0
%
% Required Global Data Structures:
%       CST
%       NRCIM
%       RES
%
% Required Data Files:
%      none.
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% static char rcsid[] = "$Id: nrSciMeritPlotSegDisp.m,v 1.2 2012/07/24 17:22:36 roberts Exp $";
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

function [xyuv h] = nrSciMeritPlotSegDisp(pmsp,plotflag)

global CST;
global NRCIM;
global RES;

M1_surfnum = 1;
NumSegs = NRCIM.SegsPerSurf(M1_surfnum);

idx = 1;
% Note that pmsp is already filtered for M1 surfaces only.  We need to find
% the CSTs that correspond to those surfaces.  This is done by filtering
% for only the M1 surfaces

xyuv = zeros(NumSegs,4);
for ii = 1:CST.Num
    if CST.Meta(ii).NrcimSurf == M1_surfnum
        % Get the transform to convert a point in local M coordinates to
        % the global M1CS.
        T = vlCsMult(vlCsInv(RES.CST_Z_TMTM1),CST.ZM(:,:,ii));
        % x and y are in Z coordinates which are taken directly from T
        xyuv(idx,1) = T(1,4);
        xyuv(idx,2) = T(2,4);
        % want to convert the dX and dY positions from M to Z
        PZ = vlCsPMult(T,[pmsp(idx,3) pmsp(idx,4) 0]');
        % to get U and V we need to subtract off the x and y
        xyuv(idx,3) = PZ(1) - xyuv(idx,1);
        xyuv(idx,4) = PZ(2) - xyuv(idx,2);
        idx = idx+1;
    end
end

if plotflag
    figure;
    h = quiver(xyuv(:,1),xyuv(:,2),xyuv(:,3),xyuv(:,4));
    title('Vector Plot of Segment Displacement (M1CS)');
else
    h = 0;
end

% end of nrSciMeritPlotSegDisp
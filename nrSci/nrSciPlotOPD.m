% File: nrSciPlotOPD.m
%
% Syntax: nrSciPlotOPD
%
% Description:
%       Plots OPD for the NRCIM version of the Merit Function Routine
%       (nrSciMerit)
%
% Input Parameters:
%       None.
%
% Output Parameters:
%       None
%
% Required Global Data Structures:
%       OC  - The OC data structure is required.
%       RES - The Merit Function result data structure.  This is
%       initialized in nrSciMerit
%
% Required Data Files:
%       N/A
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: nrSciPlotOPD.m,v 1.5 2012/10/29 22:23:15 roberts Exp $
%
% INDENT-OFF*
% $Log: 
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

function nrSciPlotOPD

global OC
global RES

if RES.PlotOPD || RES.PlotSegDisp
    fprintf('....Plotting Results\n');
else
    % zero out uncalculated results that will be reported to Ansys
    RES.PRIM_DEC_MAX = 0;
    RES.SEG_MAX_SHEAR = 0;
    RES.SEG_MAX_GAP = 0; 
    RES.SEG_MIN_GAP = 0;
end

if RES.PlotOPD 
    % plot the original uncorrected OPD
    opd = vlLmRunLom(RES.UncorrectedSp);
    rms = vlOpWtoRMS(opd);
    local_plot_opd(opd,'Unaligned telescope OPD',rms,[RES.outFileRoot '_OPD'])
    
    % plot the original uncorrected OPD of the primary only
    opd = vlLmRunLom(RES.UncorrectedSpPrimary);
    rms = vlOpWtoRMS(opd);
    local_plot_opd(opd,'Unaligned Primary Mirror OPD',rms,[RES.outFileRoot '_PrimaryOPD'])

    % plot OPDzero
    rms = vlOpWtoRMS(OC.LOM.OPDzero);
    local_plot_opd(OC.LOM.OPDzero,'Residual OPD from optical prescription [OPDzero]',rms,[RES.outFileRoot '_OPDzero'])

    % Plot the opd of the decentered and rotated primary mirror
    local_plot_opd(RES.opd_decrot,'OPD of primary mirror decenter and rotation',RES.IQ_DECROT_RMS,[RES.outFileRoot '_DecRot'])

    % Plot the opd of the decentered primary mirror
    local_plot_opd(RES.opd_dec,'OPD of primary mirror decenter',RES.IQ_DEC_RMS,[RES.outFileRoot '_Dec'])

    % Plot the opd of the rotated primary mirror
    local_plot_opd(RES.opd_rot,'OPD of primary mirror rotation',RES.IQ_ROT_RMS,[RES.outFileRoot '_Rot'])

%     % Plot the final compensated OPD
%     local_plot_opd(RES.CorrOPD,'Final Compensated OPD',RES.IQ_ALL_RMS,[RES.outFileRoot '_AllComp'])
end



if RES.PlotSegDisp
    % All M1 plots should be in the M1 coordinate system.  Therefore xyuv
    % below should be in that system.
    fprintf('\tPlotting Global Primary Mirror Segment Displacement\n');
    [xyuv h] = nrSciMeritPlotSegDisp([zeros(length(RES.sp_compensated_M1),2) RES.sp_compensated_M1],RES.PlotSegDisp);
    title(sprintf('Global M1 Segment Disp. (M1CS), Max = %.2e [m]',RES.PRIM_DEC_MAX));
    xlim([-OC.EPD/2 OC.EPD/2]);
    ylim([-OC.EPD/2 OC.EPD/2]);
    vlUtResizeFigure(gcf,700,700);
    saveas(h,[RES.outFileRoot 'SegDispGlobal']);
    saveas(h,[RES.outFileRoot 'SegDispGlobal'],'tif');
    
    fprintf('\tPlotting Relative Segment Disp\n');
    [h RES.SEGDISP_XYUV RES.SEG_MAX_SHEAR RES.SEG_MAX_GAP RES.SEG_MIN_GAP RES.SEG_GAP_SHEAR_LIST] = ...
            vlSciMeritPlotSegRelDisp(xyuv,OC.SegSize+0.3,0.4,RES.PlotSegDisp);
    title(sprintf('Relative M1 Seg. Disp. (M1CS), [Max Shear,Gap(+/-)] = %.2e,%.2e,%.2e [m]',...
            RES.SEG_MAX_SHEAR,RES.SEG_MAX_GAP,RES.SEG_MIN_GAP));
    xlim([-OC.EPD/2 OC.EPD/2]);
    ylim([-OC.EPD/2 OC.EPD/2]);
    vlUtResizeFigure(gcf,700,700);
    
    saveas(h,[RES.outFileRoot 'SegDispRelative']);
    saveas(h,[RES.outFileRoot 'SegDispRelative'],'tif');
    % Create a data structure for segment rotation as a function of
    % position
    RES.SEGROT_XYRZ = zeros(size(RES.SEGDISP_XYUV,1),3);
    RES.SEGROT_XYRZ(:,1:2) = RES.SEGDISP_XYUV(:,1:2);
    % remove zero rows from pmsp
    pmsp = [zeros(length(RES.sp_compensated_M1),2) RES.sp_compensated_M1];
    RES.SEGROT_XYRZ(:,3) = pmsp(:,8);
else
    % Calculate segment displacement results whether or not the plot is created
    % since these are reported in results in either case. The final argument to
    % these functions specifies whether the plot is created or not.
    [xyuv h] = nrSciMeritPlotSegDisp([zeros(length(RES.sp_compensated_M1),2) RES.sp_compensated_M1],RES.PlotSegDisp);
    [h RES.SEGDISP_XYUV RES.SEG_MAX_SHEAR RES.SEG_MAX_GAP RES.SEG_MIN_GAP RES.SEG_GAP_SHEAR_LIST] = ...
            vlSciMeritPlotSegRelDisp(xyuv,OC.SegSize+0.3,0.4,RES.PlotSegDisp);
end

if RES.PlotACT
    fprintf('\tPlotting Histogram of M1 Actuator Motions\n');
    h = vlSciPlotHist(RES.ActStroke,20,'Actuator Motions');
    vlUtResizeFigure(h,800,800);
    saveas(h,[RES.outFileRoot 'ActHistogram']);
    saveas(h,[RES.outFileRoot 'ActHistogram'],'tif');   
    
    % Need actuator positions in the M1CS
    M1ActPos = zeros(size(RES.ActNodePosList));
    % Calculate the transform to convert a point from Z to M1CS coordinates
    T = vlCsInv(RES.CST_Z_TMTM1);
    
    for ii=1:length(RES.ActNodePosList)
        P = vlCsPMult(T,RES.ActNodePosList(ii,:)');
        M1ActPos(ii,:) = P';
    end
        
    fprintf('\tPlotting Actuator Stroke Requirement for this load case\n');
    h = vlSciPlotActPos(M1ActPos,RES.ActStroke,'M1 Actuator Stroke in M1CS');
    vlUtResizeFigure(gcf,700,650);
    saveas(h,[RES.outFileRoot 'ActStroke']);
    saveas(h,[RES.outFileRoot 'ActStroke'],'tif');
end

fprintf('....done in %.1f Seconds\n\n',toc);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%         local_plot_opd
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function local_plot_opd(opd,titlename,rms,outfilename)
global OC

fprintf('\tPlotting %s\n',titlename);
figure;
% remove bias in opd
offset = (max(opd(OC.LOM.IAD>0))+min(opd(OC.LOM.IAD>0)))/2;
opd = opd - offset*ones(size(opd)).*OC.LOM.IAD;

opd = fliplr(flipud(reshape(opd,OC.LOM.N,OC.LOM.N)));
h = surf(opd, 'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
axis([0 OC.LOM.N 0 OC.LOM.N]);
title(sprintf('%s, RMS = %.4e',titlename,rms));
view(2);
colorbar;
vlUtResizeFigure(gcf,700,650);
saveas(h,outfilename)
saveas(h,outfilename,'tif')

return;

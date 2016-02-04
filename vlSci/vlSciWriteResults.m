% File: <vlSciWriteResults.m>
%
% Syntax: vlSciWriteResults(varargin) or vlSciWriteResults('inputfile') or vlSciWriteResults
%
% Description:
%       This function writes results from the Merit Function res structure
%       to 3 places:
%           1) Writes a verbose descripton of the results to the command
%           window.
%           2) Writes the same description to a file
%           3) Writes a text file of results for reading into Ansys
%
% Input Parameters:
%       Input arguments can either be entered on the command line or by
%       using an input file.  If using an input file, enter the arguments
%       in an ascii text file using spaces or line breaks to separate the
%       arguments.  In this case, nrSciMerit takes only one argument, the
%       name of the input text file.
%
%       If no input arguments are given, the routine will look for a file
%       called 'MatlabCmdLine.txt' and execute commands from there.
%       
%       CASEFLAG - either '1Case' or 'Gravity' 
%           '1Case' - a single load case is reported on (valid for thermal,
%           or applied load type analysis.
%           'Gravity' - analysis of the difference between gravity loading
%           at two separate elevation angles.
%       If '1Case' the BASISCASE argument is supplied
%       If 'Gravity' the PERT_NOMACT and PERT_FULL arguments are supplied
%
%       BASISCASE - the _Data.mat file containing the RES data structure
%       for the case to be analyzed.  If a 1Case analysis this is the only
%       result file supplied.  For a Gravity analysis this is the case for
%       which the optical system is assumed to be perfectly aligned.
%
%       PERT_NOMACT - Only supplied for 'Gravity' analysis cases. This is the
%       _Data.mat file for the perturbed case that is compared to the basis
%       case using the -NomAct flag.
%
%       PERT_FULL - Only supplied for 'Gravity' analysis cases.  This is
%       the _Data.mat file for the perturbed case elevation angle using
%       full gravity results.
%
%       Flags
%           '-PlotACT' - creates histogram and stroke plots
%           of required actuator motions.
%           '-PlotOPD' - creates plots of the optical path differences.
%           '-OutputSP' - outputs the SP vector to a file
%           '-PlotTriads' - plots the orientation of the
%           triads on the primary mirror.
%           '-PlotSegDisp' - plots the displacement of
%           primary mirror segments and a histogram plot of
%           the primary mirror actuator motions
%           '-MapOPD' - makes a test OPD plot to determine
%           the orientation relative to the primary.
%           '-PlotAll' - equivalent to '-PlotOPD','-PlotSegDisp','-PlotAct'
%           '-ExitMatlab' - exit Matlab after execution. This is
%                           used to exit Matlab after a call from Ansys.
%                           Default state is not to exit.

%
%
% Output Parameters:
%        None
%
% Required Global Data Structures:
%       none.
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
% $Id: vlSciWriteResults.m,v 1.18 2012/11/02 17:50:49 roberts Exp $
%
% INDENT-OFF*
% $Log: vlSciWriteResults.m,v $
% Revision 1.18  2012/11/02 17:50:49  roberts
% updated text of report messages
%
% Revision 1.17  2012/10/30 15:59:45  roberts
% Removed reporting of RES.ActFileName
%
% Revision 1.16  2012/10/29 22:26:00  roberts
% added to report that FOS CS rotates with elevation angle
%
% Revision 1.15  2012/08/03 00:48:54  roberts
% error in where report file was being written fixed
%
% Revision 1.14  2012/08/02 23:06:59  roberts
% added path to ImSource in file
%
% Revision 1.13  2012/08/02 20:25:49  roberts
% added ability to read inputs from a command file
%
% Revision 1.12  2012/07/24 20:06:16  roberts
% major update to enable reporting of both '1Case' and 'Gravity' RES structures
%
% Revision 1.11  2010/12/20 18:24:03  roberts
% removed output of untested results
%
% Revision 1.10  2010/02/22 01:39:30  roberts
% updated for MFR8.1
%
% Revision 1.9  2009/09/13 22:29:40  roberts
% removed reporting of LOS results.
%
% Revision 1.8  2008/05/14 21:27:22  roberts
% updated to report info about -NomAct flag
%
% Revision 1.7  2007/02/24 01:45:00  msmith
% Added printing of Actuator Fitting Method.
%
% Revision 1.6  2006/04/13 23:10:04  roberts
% was not correctly writing outFileRoot
%
% Revision 1.5  2006/04/11 22:26:50  roberts
% Changed primary mirror coordinates from ZZP to EXYZM1P
%
% Revision 1.4  2006/01/13 00:16:27  roberts
% cleaned up.
%
% Revision 1.3  2006/01/08 22:56:56  roberts
% numerous updates
%
% Revision 1.2  2006/01/06 20:40:42  roberts
% Modified output to Ansys
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vlSciWriteResults(varargin)

% declare integrated model global variables 
global IM
global CST
global OC
global TS
global RES
global MA       % used by other routines but loaded by this function
global NRCIM    % ditto
global SD       % ditto
global SEN
global TMTCS

tic;

if length(varargin) < 2
    if isempty(varargin)
    % read input values from MatlabCmdLine.txt file
        fid = fopen('MatlabCmdLine.txt');
        if fid == -1
            error('Unable to open the input file MatlabCmdLine.txt');
        end
    else
        % input values will be read from a file
        fid = fopen(varargin{1});
        if fid == -1
            error('Unable to open the input file %s',varargin{1});
        end
    end
    ii = 1;
    while 1 
        tmp = fscanf(fid,'%s',1);
        if isempty(tmp)
            fclose(fid);
            break;
        end
        input{ii} = tmp;
        ii = ii+1;
    end
else
    input = varargin;
end

if length(input) < 2
    error('Insufficient number of arguments supplied');
end

switch(input{1})
    case '1Case'
        fprintf('Single Load Case\n');
        RES = load(input{2});
        RES = RES.RES;
        if RES.NomAct
            error('For 1Case analysis RES.NomAct cannot be set');
        end
        ii = 3;
    case 'Gravity'
        fprintf('Gravity Analysis Case\n');
        pert_nomact_RES = load(input{2});
        pert_nomact_RES = pert_nomact_RES.RES;
        RES = pert_nomact_RES;
        if ~pert_nomact_RES.NomAct
            error('For Gravity Analysis, perturbation NomAct case RES.NomAct must be set');
        end
        pert_full_RES = load(input{3});
        pert_full_RES = pert_full_RES.RES;
        if pert_full_RES.NomAct
            error('For Gravity Analysis, perturbation full case RES.NomAct must not be set');
        end
        
        fprintf('Loading basis case, %s\n',RES.NomActPosFile);
        basisRES = load(RES.NomActPosFile);
        basisRES = basisRES.RES;
        if basisRES.NomAct
            error('For Gravity Analysis Basis case cannot have RES.NomAct set');
        end       
        
        ii = 4;
    otherwise
        error(['vlSciWriteResults: Unknown analysis flag ' input{1}]);
end

% load the NRCIM configuration file
load(RES.MatFileName);

% Set default plotting options
RES.PlotOPD = 0;
RES.PlotACT = 0;
RES.PlotTriads = 0;
RES.PlotSegDisp = 0;
RES.MapOPD = 0;
ExitMatlab = 0;

while ii <= length(input)
    switch input{ii}      
        case '-PlotAll'
            RES.PlotOPD = 1;
            RES.PlotACT = 1;
            RES.PlotSegDisp = 1;
        case '-MapOPD'
            RES.MapOPD = 1;
        case '-PlotOPD'
            RES.PlotOPD = 1;  
        case '-PlotACT'
            RES.PlotACT = 1;
        case '-OutputSP'
            RES.OutputSP = 1;   
        case '-PlotTriads' 
            RES.PlotTriads = 1;
        case '-PlotSegDisp' 
            RES.PlotSegDisp = 1;
        case '-ExitMatlab'
            ExitMatlab = 1;
        otherwise
            error(['vlSciWriteResults: Unknown input argument' input{ii}]);
    end
    % Advance to next parameter
    ii = ii + 1;    
end

% Add IMSource to the IM source file path
if isdir(fullfile('.','ImSource'))
    fprintf('\tFound IMSource sub-directory, adding to path...\n');
    addpath(fullfile(pwd,'ImSource'));
    if isdir(fullfile('.','ImSource','Data'))
        fprintf('\tFound ImSource\\Data sub-directory, adding to path...\n');
        addpath(fullfile(pwd,'ImSource','Data'));
    end
end

if RES.PlotTriads
    fprintf('....Plotting Triads...\n');
    [xyuv h] = nrCsPlotTriads(1);
    title(sprintf('Plot of Triad Orientations'));
    saveas(h,[RES.outFileRoot 'Triads']);
    saveas(h,[RES.outFileRoot 'Triads'],'tif');
end

if RES.MapOPD  
    % test orientation of the OPD plot
    sp = zeros(CST.Num,6);

    segspercon = NRCIM.SegsPerSurf(1)/6;
    for con = 1:6
      for seg = 1:segspercon
          sp((con-1)*segspercon+seg,3) = -((con-1)*segspercon+seg)*1e-6; % piston
      end
    end

    sp = sp';
    opd = vlLmRunLom(sp(:));
    % reshape to square array and flip around both axes
    opd = rot90(reshape(opd,OC.LOM.N,OC.LOM.N),2);
    figure;
    h = surf(opd, 'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
    axis([0 OC.LOM.N 0 OC.LOM.N]);
    title(sprintf('Test of OPD orientation'));
    view(2);
    colorbar;
    saveas(h,[RES.outFileRoot 'OPD_Orient']);
    saveas(h,[RES.outFileRoot 'OPD_Orient'],'tif');
end




%% create plots
nrSciPlotOPD;

%% Open output files
% open the output file for the Merit Function Run report
report_fid = fopen(RES.ReportFileName,'wt');
if report_fid == -1
    error('Unable to open Report file %s\n',RES.ReportFileName);
end

% open the output file for the Ansys result file
result_fid = fopen(RES.ResultFileName,'wt');
if result_fid == -1
    error('Unable to open Result file %s\n',RES.ResultFileName);
end

%% Write results
% write metadata report to the screen
localWriteMetaData(RES,1); 
% write metadata report to file 
localWriteMetaData(RES,report_fid);
% write Ansys result to file
localWriteAnsysResult(RES,result_fid);

% write actuator motion report to the screen
localWriteActResults(RES,1); 
% write actuator motion report to file 
localWriteActResults(RES,report_fid);

if RES.NomAct % Gravity Case
    if isempty(basisRES.ElevAng)
        error('Basis case RES.ElevAng is empty');
    end
    if isempty(RES.ElevAng)
        error('RES.ElevAng is empty');
    end
    GRAVRES = nrSciCalcMotion2Cases(basisRES,pert_full_RES,pert_full_RES.ElevAng-basisRES.ElevAng,1);
    nrSciCalcMotion2Cases(basisRES,pert_full_RES,pert_full_RES.ElevAng-basisRES.ElevAng,report_fid);
    
    tempRES = RES;
    RES = pert_nomact_RES;
    save(strcat(RES.outFileRoot,'_Data'),'RES','GRAVRES');
    RES = tempRES;
else
    localWriteReport1Case(RES,1); 
    % write report to file 
    localWriteReport1Case(RES,report_fid);
end

if RES.OutputSP
    sp = RES.UncorrectedSp;
	save([RES.outFileRoot 'SP.txt'],'sp','-ASCII');
end

% close the files
status = fclose(result_fid);
if status == -1
    error( [ 'Could not close file ' RES.ResultFileName '.' ] );
end

status = fclose(report_fid);
if status == -1
    error( [ 'Could not close file ' RES.ReportFileName '.' ] );
end

% Check for ExitMatlab flag and exit if true
if ExitMatlab 
    exit;
else
    return;
end

% end of vlSciWriteResults
end

function localWriteMetaData(RES,fid)
fprintf(fid,'\n\n\n');
fprintf(fid,'************************************************************************\n');
fprintf(fid,'******           Herzberg Institute of Astrophysics               ******\n');
fprintf(fid,'******      Astronomy Technology Research Group - Victoria        ******\n');
fprintf(fid,'******                                                            ******\n');
fprintf(fid,'******        Telescope Structure Merit Function Tool             ******\n');
fprintf(fid,'******                                                            ******\n');
fprintf(fid,'************************************************************************\n');
fprintf(fid,'* (c) <2005-2012>				        (c) <2005-2012>                 *\n');
fprintf(fid,'* National Research Council         Conseil national de recherches     *\n');
fprintf(fid,'* Ottawa, Canada, K1A 0R6 		    Ottawa, Canada, K1A 0R6            *\n');
fprintf(fid,'* All rights reserved			    Tous droits reserves               *\n');
fprintf(fid,'************************************************************************\n');

fprintf(fid,'\nGeneral Run Data...\n');
fprintf(fid,'\t...Input and Output Files\n');
fprintf(fid,'\t\tInput Displacement File Name [YFileName] is: %s\n',RES.YFileName);
fprintf(fid,'\t\tInput Matlab Data Structure File Name [MatFileName] is: %s\n',RES.MatFileName);
fprintf(fid,'\t\tRoot File Name for Output Files [outFileRoot] is: %s\n',RES.outFileRoot);
fprintf(fid,'\t\tActuator Fitting Method [ActFit] is: %s\n',RES.ActFit);
fprintf(fid,'\t\tOutput Ansys Result File Name [ResultFileName] is: %s\n',RES.ResultFileName);

fprintf(fid,'\n\t...Units...\n');
fprintf(fid,'\t\tUnless otherwise stated all results are in meters and degrees\n');

fprintf(fid,'\n\t...General Information\n');
fprintf(fid,'\t\tFile Directory [Dir] is: %s\n',RES.DIR);
fprintf(fid,'\t\tDate and Time of Run [Time] is: %s\n',RES.Time);
fprintf(fid,'\t\tRun Execution Time [ExecTime] is: %.1f [Seconds]\n',RES.ExecTime);
fprintf(fid,'\t\tMatlab Version [Ver] is: %s\n',RES.Ver);

fprintf(fid,'\n\t...Actuator Stroke Calculation Data\n');
fprintf(fid,'\t\tMethod of Actuator Stroke Minimization [ActFit] is: %s\n',RES.ActFit);
if RES.NomAct
    fprintf(fid,'\t\tActuator Stroke is Calculated Relative to another load case\n');
    fprintf(fid,'\t\tBasis Load Case File [NomActPosFile] is : %s\n',RES.NomActPosFile);
else
    fprintf(fid,'\t\tActuator Stroke is Calculated Relative to no load case (zero thermal & gravity)\n');  
end
fprintf(fid,'\n');

return
% end of localWriteMetaData
end

function localWriteActResults(RES,fid)
fprintf(fid,'Actuator Stroke and Primary Mirror Segment Information\n');

fprintf(fid,'\t...Required Primary Mirror Actuator Motions\n');
fprintf(fid,'\t\tMaximum Required Actuator Stroke [ACT_MAX] is: %.3e\n',RES.ACT_MAX);
fprintf(fid,'\t\tMean Actuator Stroke [ACT_MEAN] is: %.3e\n',RES.ACT_MEAN);
fprintf(fid,'\t\tRMS Actuator Stroke [ACT_RMS] is: %.3e\n',RES.ACT_RMS);
fprintf(fid,'\n');
    
fprintf(fid,'Primary Mirror Residual Segment Motions after alignment...\n');
fprintf(fid,'\t...Clocking\n');
fprintf(fid,'\t\tMaximum segment clocking [PRIM_ROT_MAX] is %.3e[deg]\n',RES.PRIM_ROT_MAX);
fprintf(fid,'\t\tMean segment clocking [PRIM_ROT_MEAN] is %.3e [deg]\n',RES.PRIM_ROT_MEAN);
fprintf(fid,'\t\tRMS segment clocking [PRIM_ROT_RMS] is %.3e [deg]\n',RES.PRIM_ROT_RMS);

fprintf(fid,'\t...Decenter magnitude\n');
fprintf(fid,'\t\tMaximum segment decenter [PRIM_DEC_MAX] is %.3e[m]\n',RES.PRIM_DEC_MAX);
fprintf(fid,'\t\tMean segment decenter [PRIM_DEC_MEAN] is %.3e [m]\n',RES.PRIM_DEC_MEAN);
fprintf(fid,'\t\tRMS segment decenter [PRIM_DEC_RMS] is %.3e [m]\n',RES.PRIM_DEC_RMS);

if RES.PlotSegDisp
    fprintf(fid,'\t...Relative Primary Mirror Segment Motion\n');  
    fprintf(fid,'\t\tMaximum relative shear between segment edges [SEG_MAX_SHEAR] is %.3e[m]\n',RES.SEG_MAX_SHEAR);
    fprintf(fid,'\t\tLargest increase in relative segment gap distance [SEG_MAX_GAP] is %.3e [m]\n',RES.SEG_MAX_GAP);
    fprintf(fid,'\t\tLargest decrease in relative segment gap distance [SEG_MIN_GAP] is %.3e [m]\n',RES.SEG_MIN_GAP);
end

fprintf(fid,'\nImage quality results due to segment motions after alignment (without warping harness adjustments)\n');
fprintf(fid,'\t...Image quality degradation due to residual decenter (Dx,Dy) and clocking (Rz) of primary mirror segments\n');
fprintf(fid,'\t\tRMS Wavefront [IQ_DECROT_RMS] is: %.3e [m]\n',RES.IQ_DECROT_RMS);

fprintf(fid,'\t...Image quality degradation due to residual decenter (Dx,Dy) of the primary mirror segments\n');
fprintf(fid,'\t\tRMS Wavefront [IQ_DEC_RMS] is: %.3e [m]\n',RES.IQ_DEC_RMS);

fprintf(fid,'\t...Image quality degradation due to residual clocking (Rz) of the primary mirror segments\n');
fprintf(fid,'\t\tRMS Wavefront [IQ_ROT_RMS] is: %.3e [m]\n',RES.IQ_ROT_RMS);
end

function localWriteReport1Case(RES,fid)
fprintf(fid,'\nTelescope Alignment Information\n');

fprintf(fid,'\t...Primary Alignment Information\n');
fprintf(fid,'\t\tCenter of Curvature of aligned primary mirror relative to original coordinates\n');
fprintf(fid,'\t\t\tX,Y,Z Coordinates [COC_X,COC_Y,COC_Z] are: %.3e,%.3e,%.6e [m]\n',RES.COC_X,RES.COC_Y,RES.COC_Z);
fprintf(fid,'\n');

fprintf(fid,'\t...Primary Mirror Pointing after Alignment and Phasing\n');
fprintf(fid,'\t\tPrimary mirror orientation relative to original coordinates [EXYZM1P]\n');
fprintf(fid,'\t\t\tEuler angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',RES.EXYZM1P);
fprintf(fid,'\n');
fprintf(fid,'\t...Uncorrected Secondary Mirror Displacement in Secondary Mirror Coordinates\n');
fprintf(fid,'\t\tSecondary mirror position from M2 alignment coordinates [M2A_M2D]\n');
fprintf(fid,'\t\t\tEuler angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',RES.EULM2A_M2D);
fprintf(fid,'\t\tSecondary mirror position from M2 original coordinates [M2O_M2D]\n');
fprintf(fid,'\t\t\tEuler angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',RES.EULM2O_M2D);
fprintf(fid,'\n');
fprintf(fid,'\t...Uncorrected Secondary Mirror Position in Primary Mirror Coordinates\n');
fprintf(fid,'\t\tSecondary mirror position from M1 alignment coordinates [M1A_M2D]\n');
fprintf(fid,'\t\t\tEuler angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',RES.EULM1A_M2D);
fprintf(fid,'\t\tSecondary mirror position from M1 original coordinates [M1O_M2D]\n');
fprintf(fid,'\t\t\tEuler angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',RES.EULM1O_M2D);
fprintf(fid,'\n');
fprintf(fid,'\t...Uncorrected Tertiary Mirror Displacement in Tertiary Mirror Coordinates\n');
fprintf(fid,'\t\tTertiary mirror position from  M3 alignment coordinates [M3P]\n');
fprintf(fid,'\t\t\tEuler angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',RES.EULM3A_M3D);
fprintf(fid,'\t\tTertiary Mirror position from  M3 original coordinates [M3N]\n');
fprintf(fid,'\t\t\tEuler angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',RES.EULM3O_M3D);
fprintf(fid,'\n');
fprintf(fid,'\t...Uncorrected Tertiary Mirror Position in Primary Mirror Coordinates\n');
fprintf(fid,'\t\tTertiary mirror position from M1 alignment coordinates [M3P]\n');
fprintf(fid,'\t\t\tEuler angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',RES.EULM1A_M3D);
fprintf(fid,'\t\tTertiary Mirror position from M1 original coordinates [M3N]\n');
fprintf(fid,'\t\t\tEuler angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',RES.EULM1O_M3D);
fprintf(fid,'\n');
fprintf(fid,'\t...MFR Focal Plane Coordinate System Rotates with Elevation Angle\n');
fprintf(fid,'\t...Uncorrected Focal Plane Displacement in Focal Plane Coordinates\n');
fprintf(fid,'\t\tDisplaced Focal Plane position from FOC alignment coordinates[FOSP]\n');
fprintf(fid,'\t\t\tEuler angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',RES.EULFOSA_FOSD);
fprintf(fid,'\t\tDisplaced Focal Plane position from FOC original coordinates [FOSN]\n');
fprintf(fid,'\t\t\tEuler angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',RES.EULFOSO_FOSD);
fprintf(fid,'\n');
fprintf(fid,'\t...Uncorrected Focal Plane Position in Primary Mirror Coordinates\n');
fprintf(fid,'\t\tDisplaced Focal Plane position from M1 alignment coordinates[FOSP]\n');
fprintf(fid,'\t\t\tEuler angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',RES.EULM1A_FOSD);
fprintf(fid,'\t\tDisplaced Focal Plane position from M1 original coordinates [FOSN]\n');
fprintf(fid,'\t\t\tEuler angles (DX,DY,DZ,RX,RY,RZ) are...\n');
fprintf(fid,'\t\t\t\t%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',RES.EULM1O_FOSD);
fprintf(fid,'\n');
% fprintf(fid,'\t...Required system compensation to center and focus image\n');
% fprintf(fid,'\t\tTertiary Mirror Tip/Tilt [TER_POINT_RX,RY] = %.3e,%.3e [deg]\n',...
%                 RES.TER_POINT_RX,RES.TER_POINT_RY);
% fprintf(fid,'\t\tSecondary Mirror Focus [SEC_FOCUS] = %.3e [m]\n',RES.SEC_FOCUS);


% fprintf(fid,'\t...Residual image quality through entire system after refocus and LOS correction\n');
% fprintf(fid,'\t\tRMS Wavefront [IQ_ALL_RMS] is: %.3e [m]\n',RES.IQ_ALL_RMS);
% fprintf(fid,'\t\t50%% Encircled Energy [IQ_ALL_50EE] is: %.3e [arcsec]\n',RES.IQ_ALL_50EE);
% fprintf(fid,'\t\t80%% Encircled Energy [IQ_ALL_80EE] is: %.3e [arcsec]\n',RES.IQ_ALL_80EE);
% fprintf(fid,'\t\tStrehl Ratio [IQ_ALL_STREHL] is: %.3e\n',RES.IQ_ALL_STREHL);

% fprintf(fid,'\t...Line of Sight Error through System Before Correction\n');
% fprintf(fid,'\t\tLOS position at focal plane X,Y [LOS_X,Y] is: %.3e,%3e [m]\n',RES.LOS_X,RES.LOS_Y);
% 
% fprintf(fid,'\t...Plate Scale\n');
% fprintf(fid,'\t\tPlate Scale of nominal optical design is %.6e [m/arcsec]\n',RES.PlateIdeal);
% fprintf(fid,'\t\tPlate Scale of as-aligned telescope due to refocus is %.6e [m/arcsec]\n',RES.PlateAlign);
% fprintf(fid,'\t\tPlate Scale difference is %.3e [m/arcsec]\n',RES.PlateIdeal-RES.PlateAlign);

% Zernikes commented out until verified
% fprintf(fid,'\t...Zernike Terms for complete system after LOS and Focus Correction\n');
% fprintf(fid,'\t\tZER11C:\t%+.4e\n',RES.ZER11C);
% fprintf(fid,'\t\tZER11S:\t%+.4e\n',RES.ZER11S);
% fprintf(fid,'\t\tZER20:\t%+.4e\n',RES.ZER20);
% fprintf(fid,'\t\tZER22C:\t%+.4e\n',RES.ZER22C);
% fprintf(fid,'\t\tZER22s:\t%+.4e\n',RES.ZER22S);
% fprintf(fid,'\t\tZER31C:\t%+.4e\n',RES.ZER31C);
% fprintf(fid,'\t\tZER31S:\t%+.4e\n',RES.ZER31S);
% fprintf(fid,'\t\tZER33C:\t%+.4e\n',RES.ZER33C);
% fprintf(fid,'\t\tZER33S:\t%+.4e\n',RES.ZER33S);
% fprintf(fid,'\t\tZER40:\t%+.4e\n',RES.ZER40);
% fprintf(fid,'\t\tZER42C:\t%+.4e\n',RES.ZER42C);
% fprintf(fid,'\t\tZER42S:\t%+.4e\n',RES.ZER42S);
% fprintf(fid,'\t\tZER44C:\t%+.4e\n',RES.ZER44C);
% fprintf(fid,'\t\tZER44S:\t%+.4e\n',RES.ZER44S);

fprintf(fid,'\n\n\n');

return
% end of localWriteReport
end

function localWriteAnsysResult(RES,fid)

local_printline(fid,6,RES.ACT_MAX,RES.ACT_RMS,RES.ACT_MEAN);
local_printline(fid,6,RES.SEG_MAX_SHEAR,RES.SEG_MAX_GAP,RES.SEG_MIN_GAP);
local_printline(fid,6,RES.PRIM_ROT_MAX,RES.PRIM_ROT_MEAN,RES.PRIM_ROT_RMS);
local_printline(fid,6,RES.PRIM_DEC_MAX,RES.PRIM_DEC_MEAN,RES.PRIM_DEC_RMS);
%local_printline(fid,6,RES.LOS_X,RES.LOS_Y);
%local_printline(fid,6,RES.IQ_ALL_RMS,RES.IQ_DECROT_RMS,RES.IQ_DEC_RMS,RES.IQ_ROT_RMS);
%local_printline(fid,6,RES.EULZZP);
%local_printline(fid,6,RES.EULM2P);
%local_printline(fid,6,RES.EULM2N);
%local_printline(fid,6,RES.EULM3P);
%local_printline(fid,6,RES.EULM3N);
%local_printline(fid,6,RES.EULFOSP);
%local_printline(fid,6,RES.EULFOSN);

return;
% end of localWriteAnsysResult
end

function local_printline(fid,NumPerLine,varargin)

NumPrinted = 0;
NumInputs = nargin - 2;

for i = 1:NumInputs
    for j = 1:length(varargin{i})
        fprintf(fid,'%17.8e',varargin{i}(j));
    end
    NumPrinted = NumPrinted + j;
end

for i = NumPrinted+1:NumPerLine
    fprintf(fid,'%17.8e',0.0);
end

fprintf(fid,'\n');

% end of local_printline
end

% File: nrSciMeritGetImSource.m
%
% Syntax: nrSciMeritGetImSource
%
% Description:
%       Copies Science Merit Function source files to a set of directories to allow
%       distribution of the files.  Must be run from the directory you want
%       the 'IMSourceFiles' directory to appear below, and there must be a
%       'data' subdirectory.  Path must be set to make all files
%       accessible.
%
% Input Parameters:
%       none.
%
% Output Parameters:
%       none.
%
% Required Global Data Structures:
%       none.
%
% Required Data Files:
%       none.
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: nrSciMeritGetImSource.m,v 1.5 2012/10/05 21:53:40 roberts Exp roberts $
%
% INDENT-OFF*
% $Log: nrSciMeritGetImSource.m,v $
% Revision 1.5  2012/10/05 21:53:40  roberts
% Updated for MFR10
%
% Revision 1.4  2011/01/19 00:40:43  roberts
% Included M1 segment CS output .m file
%
% INDENT-ON*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%           Herzberg Institute of Astrophysics                  %%%%%
%%%%%%      Astronomy Technology Research Group - Victoria           %%%%%
%
% (c) <2010>				        (c) <2010>
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

function nrSciMeritGetImSource

%% Check if the current directory is named 'M_Files'
CurrentDir = pwd;
if ~strcmp(CurrentDir(end-length('M_Files')+1:end),'M_Files')
    fprintf('This routine creates files relative to the current directory\n');
    fprintf('...The current directory path is \n\t\t%s\n',pwd);
    reply = input('Do you want to continue to create files relative to this path? Y/N [Y]: ','s');
    if upper(reply) ~= 'Y'
        return
    end
end


fprintf('This routine will copy the Merit Function m-files\n');
reply = input('Have you set your path to the appropriate directory tree? Y/N [Y]: ','s');
if upper(reply) ~= 'Y'
    return
end

%% nrCim functions
% Return to Top Level directory
change_directory('nrCim');

local_getfile('nrCimCalcPerturbations.m');
local_getfile('nrCimCheckFEAInterface.m');
local_getfile('nrCimCP.m');
local_getfile('nrCimCreateFeaDataScript.m');
local_getfile('nrCimND.m');

%% vlSci Files
cd('..');  
change_directory('vlSci');

local_getfile('vlSciActStroke.m');
local_getfile('vlSciMeritGetVersion.m');
local_getfile('vlSciMeritPlotSegRelDisp.m');
local_getfile('vlSciPlotActPos.m');
local_getfile('vlSciPlotActPosVel.m');
local_getfile('vlSciPlotHist.m');
local_getfile('vlSciWriteResults.m');

%% nrSci Files
cd('..');  
change_directory('nrSci');

local_getfile('nrSciMerit.m');
local_getfile('nrSciGravity.m');

local_getfile('nrSciMeritGetImSource.m');

local_getfile('nrSciAlignTelescope.m');
local_getfile('nrSciCalcMotion2Cases.m');
local_getfile('nrSciCalcResults.m');
local_getfile('nrSciGetTMTM1CSPerturbations.m');
local_getfile('nrSciMeritPlotSegDisp.m');
local_getfile('nrSciPlotOPD.m');
%% ansys functions
cd('..');  % Return to Top Level directory
change_directory('ansys');

local_getfile('vlAnFormatInt.m');           % pcode files that are in script directories

%% vlCs functions
cd('..');  % Return to Top Level directory
change_directory('vlCs');
local_getfile('vlCsCheckTrans.m');
local_getfile('vlCsEulerXYZ.m');
local_getfile('vlCsInit.m');
local_getfile('vlCsInv.m');
local_getfile('vlCsInvEulerXYZ.m');
local_getfile('vlCsMult.m');
local_getfile('vlCsPMult.m');
local_getfile('vlCsRigidFrameTrans.m');
local_getfile('vlCsRotX.m');
local_getfile('vlCsRotY.m');
local_getfile('vlCsRotZ.m');
local_getfile('vlCsTrans.m');
local_getfile('vlCsTriadTrans.m');
local_getfile('vlCsUnit.m');

%% nrCs functions
cd('..');  % Return to Top Level directory
change_directory('nrCs');
local_getfile('nrCsPlotTriads.m');

%% vlLom functions
cd('..');  % Return to Top Level directory
change_directory('vlLom');

local_getfile('vlLmRunLom.m');

%% vlOptical functions
cd('..');  % Return to Top Level directory
change_directory('vlOptical');

local_getfile('vlOpDecComp.m');
local_getfile('vlOpGetZer.m');
local_getfile('vlOpGmm.m');
local_getfile('vlOpInit.m');
local_getfile('vlOpMovetoSphere.m');
local_getfile('vlOpPSFtoEE.m');
local_getfile('vlOpPSFtoStrehl.m');
local_getfile('vlOpSag.m');
local_getfile('vlOpScalarProd.m');
local_getfile('vlOpWtoPSFFFT.m');
local_getfile('vlOpWtoRMS.m');
local_getfile('vlOpWtoZer.m');
local_getfile('vlOpZerInit.m');
local_getfile('vlOpZerNumero.m');

%% vlRaysIn functions
cd('..');  % Return to Top Level directory
change_directory('vlRaysIn');

local_getfile('vlRyLomGrid.m');

%% vlUtil functions
cd('..');  % Return to Top Level directory
change_directory('vlUtil');

local_getfile('vlUtFindFiles.m');
local_getfile('vlUtGetIntFromString.m');
local_getfile('vlUtResizeFigure.m');

%% Ansys Macros

cd('..');  % Return to Top Level directory
change_directory('AnsysMacro');

local_getfile('vlAnMeritInit.mac');
local_getfile('vlAnMeritStaticDisp.mac');
local_getfile('vlAnMeritReadResults.mac');


%% Copy model specific CST functions into Data subdirectory.
% Data functions % MFR R6 and later

cd('..');  % Return to Top Level directory
change_directory('Data');

local_getfile('TMTRCV1.m');
local_getfile('vlCsRcCSYS0toM1CS_3mirr.m');
local_getfile('vlCsRcGetDistances_3mirr.m');
local_getfile('vlCsRcZFocus_3mirr.m');
local_getfile('vlCsRcZsecondary_3mirr.m');
local_getfile('vlCsRcZtertiary_3mirr.m');

% Change directory back to the original directory
cd(CurrentDir);

return;
% end of nrSciMeritGetImSource

%%   local_getfile(fname)

function local_getfile(fname)

fileToCopy = which(fname);

if isempty(fileToCopy)
    fprintf('>>> Warning: %s not found\n',fname);
    return;
end

fileHere = dir(['.\' fname]);

if isempty(fileHere)
    fprintf('Copying %s\n',fileToCopy);
    copyfile(fileToCopy,'.');
else
    fprintf('Found %s in directory: Not Copying.\n',fname);
end


%%  change_directory(dirname)
 
function change_directory(dirname)

persistent allow_create;

% If the directory does not exist, then check with user if it should be
% created or not.

if isempty(dir(dirname))

    % If the user has not allowed us to always create the directory
    % then ask.

    if isempty(allow_create)

        % Keep asking user until they provide an acceptable answer.

        while (1)

            fprintf('The directory "%s" does not exist', dirname);

            reply = input('Create directory (Yes, No, Always)? Y/N/A [N]: ','s');

            % If no answer then supply the default.

            if isempty(reply)
                reply = 'N';
            end

            reply = upper(reply);

            if ismember(reply, [ 'Y', 'N', 'A' ])
                break;
            end
        end  % while


        % Check if the user does not want the directory created.

        if upper(reply) == 'N'
            error('Installation stopped by user');
        end


        % Set allow_create if directory creation is always allowed.

        if upper(reply) == 'A'
            allow_create = 1;
        end

    end

    mkdir(dirname);
end

% The directory was either created or it previously existed; cd into it.

cd(dirname);

% end of local_getfile

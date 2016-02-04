%   File: <vlOpInit.m> 
%
%   Syntax: vlOpInit
%
%   Description:
%       Fills OC global data structure.  All changes to the OC data
%       structure should be made here.
%       
%   Input Parameters:
%        none
%
%   Output Parameters:
%        Modifies global data structure OC
%
%   Required Global Data Structures:
%        OC, IM
%
%   Required Files:
%        none
%

%
% Revision History
%
% static char rcsid[] = "$Id: vlOpInit.m,v 1.19 2005/01/14 18:18:14 stretchn Exp $";
% INDENT-OFF*
% $Log: vlOpInit.m,v $
% Revision 1.19  2005/01/14 18:18:14  stretchn
% Changed more 'string' commands to 'num2str'
%
% Revision 1.18  2004/11/30 00:25:14  dunn
% Changed call to vlLmInitLOM to vlLmInitLom - important on Linux.
%
% Revision 1.17  2004/08/13 17:36:16  mckenzie
% Changed vlMzTimeOut from 120 to 300 seconds.
%
% Revision 1.16  2003/12/23 23:50:49  vlotim
% Fixed error message
%
% Revision 1.15  2003/12/19 01:11:26  stretchn
% Now performs LOM initialization if you are using zemax, but want to
% generate outputs that require lom values
%
% Revision 1.14  2003/12/18 23:43:44  stretchn
% Now init Zer and DM for Zemax OR LOM - because these can now be calculated with Zemax
%
% Revision 1.13  2003/12/05 23:59:53  stretchn
% Now calls vlOpDMInit - only calls this and/or vlOpZerInit when they
% are required (DM and/or zernikes are included in IM.SimOpdOuts)
%
% Revision 1.12  2003/11/24 22:09:25  stretchn
% Changed string to num2str
%
% Revision 1.11  2003/11/22 00:02:02  stretchn
% Fixed minor bug
%
% Revision 1.10  2003/11/18 18:24:10  stretchn
% Now allows .m to be included in OCFileSTr
%
% Revision 1.9  2003/10/30 21:13:28  stretchn
% Fixed loading of LomModelFile, so rest of OC isn't deleted
% added 'tic' commands so 'tocs' aren't unmatched.
%
% Revision 1.8  2003/10/24 18:36:17  stretchn
% Improved efficiency - removed unnecessary opening and closing of Zemax link.
%
% Revision 1.7  2003/09/19 18:47:29  stukasa
% Header update.
%
% Revision 1.6  2003/08/22 19:25:46  roberts
% Removed need for Zemax to be running when loading LOM
%
% Revision 1.5  2003/08/08 20:29:40  roberts
% added support for IM.LomModeFile
%
% Revision 1.4  2003/07/15 19:10:16  roberts
% Added longer timeout for linear optics model.
%
% Revision 1.3  2003/07/14 23:53:50  roberts
% Added timing functions for LOM and ZER initialization
%
% Revision 1.2  2003/07/04 17:57:01  dunn
% try/catch around zemax link.
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

function vlOpInit()

global OC
global IM

% check if the data structure exists.  If not, fail with error
try 
    filename = IM.OCFileStr;
catch
    error('Global Data Element IM.OCFileStr is not defined');
end

[path,name,ext]=fileparts(filename);

% If the path is valid, run the routine, else fail with error
try
    feval(name);
catch
    mess=sprintf('Error loading IM.OCFileStr.  Probable cause:\nOptical Configuration File %s Not Found on Path.\nError message:\n%s',IM.OCFileStr, lasterr);
    error(mess);
end

% Perform the following if the Linear Optical Model is being used.
if strcmp(IM.OpticalEngine,'All') || strcmp(IM.OpticalEngine,'LOM') || any([strcmp(IM.SimOpdOuts,'zernikes'),strcmp(IM.SimOpdOuts,'zerrms'),strcmp(IM.SimOpdOuts,'DM')])
    if ~isfield(IM,'LomModelFile') || ~exist(which(IM.LomModelFile)) % don't insist on zemax running if LOM is read from file.
        try
            vlMzOpen;
            [ZemFile, ZemPath, ZemVer] = vlMzGetZemaxFilePathVer;
            
            disp(sprintf('Zemax\n   Version = [%s]\n   Zemax Filename = [%s]\n   Zemax Path = [%s]', ...
                ZemVer,ZemFile,ZemPath));
        catch
            vlMzClose;
            disp('Failed to open the matlab to zemax link');
            rethrow(lasterror);
        end
        
        if ~strcmp(upper(ZemFile),upper(IM.ZemaxPrescription))
            vlMzClose;
            cd(currentDir);
            error(sprintf(['   IM.ZemaxPrescription = [%s]\n    Actual prescription = [%s]\n'...
                    'Zemax Prescription Mismatch'],IM.ZemaxPrescription,ZemFile));
        end     
        
        disp(['    Initializing LOM Linear optics model data structure ...']);
        % set timeout to allow large traces
        tic;
        vlMzSetTimeout(300);
        vlLmInitLom;
        vlMzClose;
        disp(['        ...Done, elapsed time = ' num2str(toc)]);
    else
        disp(['    Reading LOM Linear optics model data from file ...']);
        tempLOM = load(IM.LomModelFile,'OC');
        OC.LOM = tempLOM.OC.LOM;
        clear tempLOM;
    end
end

%Only intialize Zernikes structure if Zernikes are being output
if any([strcmp(IM.SimOpdOuts,'zernikes'),strcmp(IM.SimOpdOuts,'zerrms')])
    disp(['    Initializing ZER Zernike data structure ...']);
    tic;
    vlOpZerInit;
    disp(['        ...Done, elapsed time = ' num2str(toc)]);
end

%Only initialize DM structure if the DM representation of the OPD is
%being output
if any(strcmp(IM.SimOpdOuts,'DM'))
    disp(['    Initializing DM Deformable Mirror data structure ...']);
    tic;
    vlOpDMInit;
    disp(['        ...Done, elapsed time = ' num2str(toc)]);
end


% End of vlOpInit.m

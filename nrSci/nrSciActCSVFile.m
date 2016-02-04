% File: <nrSciActCSVFile.m>
%
% Syntax: nrSciActCSVFile(databasefile, filenamespec, actfilename, pointfilename, decfilename)
%
% Description:
%        Generates Excel .csv files that contain data to calculate
%        actuator rates, M1 Alignment coordinate systems and segment
%        decenter information as a function of elevation angle.
% 
% Input Parameters:
%       databasefile - [filename] filename of NRCIM database for the MFR
%       release.
%       filenamespec - [filename] specifies pattern to match all filenames.
%           For example, '*Align0_Data.mat'.
%       actfilename - [filename] specifies output .csv file for actuator
%       data.
%       pointfilename - [filename] specifies output .csv file for EulerXYZ
%       pointing of primary mirror at each elevation angle.
%       decfilename - [filename] specifies output .csv file for segment
%       decenters.
%
% Output Parameters:
%       Creates three .csv files, as described above
%
% Required Global Data Structures:
%         CST
%         NRCIM
%         RES
%         TMTCS
%
% Required Data Files:
%       Output Data files from MFR runs of gravity loading analysis 
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: nrSciActCSVFile.m,v 1.2 2012/11/01 22:24:32 roberts Exp $
%
% INDENT-OFF*
% $Log: nrSciActCSVFile.m,v $
% Revision 1.2  2012/11/01 22:24:32  roberts
% fixed header comments, no change to code
%
% Revision 1.1  2012/11/01 21:44:16  roberts
% Initial revision
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


function nrSciActCSVFile(databasefile,filenamespec,actfilename,pointfilename,decfilename)

%% Setup
global CST;
global NRCIM;
global RES;
global TMTCS;

% Add IMSource to the IM source file path
if isdir(fullfile('.','ImSource'))
    fprintf('\tFound IMSource sub-directory, adding to path...\n');
    addpath(fullfile(pwd,'ImSource'));
    if isdir(fullfile('.','ImSource','Data'))
        fprintf('\tFound ImSource\\Data sub-directory, adding to path...\n');
        addpath(fullfile(pwd,'ImSource','Data'));
    end
end

% load NRCIM database file
load(databasefile);

D=dir(filenamespec);
numcases = length(D);

% load first file 
load(D(1).name);

M1_surfnum = 1;
NumSegs = NRCIM.SegsPerSurf(M1_surfnum);

% assign data sizes
numacts = length(RES.ActStroke);
Elevangs = zeros(numcases,1);
Actstrokes = zeros(numacts,numcases);
Decenters = zeros(NumSegs,numcases*2);

Pointings = zeros (7,numcases); % elevation angle and EulerXYZ for each case

%% Get actuator stroke and M1 Alignment

for ii=1:length(D)
    fprintf('Found file: %s\n',D(ii).name);
    load(D(ii).name);
    
    % Elevation angles and actuator strokes
    Elevangs(ii,1) = RES.ElevAng;
    Actstrokes(:,ii) = RES.ActStroke;
    Pointings(1,ii) = RES.ElevAng;
    Pointings(2:7,ii) = vlCsEulerXYZ(GRAVRES.T_M1OP_M1AP)';
    
    % Decenters
    idx = 1;
    % Note that pmsp is already filtered for M1 surfaces only.  We need to find
    % the CSTs that correspond to those surfaces.  This is done by filtering
    % for only the M1 surfaces

    pmsp = RES.sp_compensated_M1;
    xyuv = zeros(NumSegs,4);
    M1xy = zeros(NumSegs,2);
    segnames = cell(NumSegs,1);

    for jj = 1:CST.Num
        if CST.Meta(jj).NrcimSurf == M1_surfnum
            % Get the transform to convert a point in local M coordinates to
            % the global M1CS.
            T = vlCsMult(vlCsInv(RES.CST_Z_TMTM1),CST.ZM(:,:,jj));
            % x and y are in Z coordinates which are taken directly from T
            xyuv(idx,1) = T(1,4);
            xyuv(idx,2) = T(2,4);
            M1xy(idx,1) = TMTCS.M1Seg(1,4,idx);
            M1xy(idx,2) = TMTCS.M1Seg(2,4,idx);
            % check distance between segmentation database and M CS segment
            % position
            segposerror = norm([TMTCS.M1Seg(1,4,idx)-T(1,4) TMTCS.M1Seg(2,4,idx)-T(2,4)]);
            if segposerror > 0.1
                fprintf('\nSegment position exceeds tolerance for segment %d. Offset is %f\n',idx,segposerror);
            end
            % want to convert the dX and dY positions from M to M1CS
            PZ = vlCsPMult(T,[pmsp(idx,1) pmsp(idx,2) 0]');
            % to get U and V we need to subtract off the x and y
            xyuv(idx,3) = PZ(1) - xyuv(idx,1);
            xyuv(idx,4) = PZ(2) - xyuv(idx,2);
            segnames{idx,1} = CST.Meta(jj).ANameSurf;
            idx = idx+1;
        end
    end
    col = 1+(ii-1)*2;
    Decenters(:,col:col+1) = xyuv(:,3:4);
end

%% Reorder data by elevation angle
[Elevangs,Order] = sort(Elevangs);

Pointings = Pointings(:,Order');
Actstrokes = Actstrokes(:,Order');

% Generate vector to reorder decenters
deco = 1+(Order-1)*2;
deco=[deco deco+1]';
deco=deco(:);
Decenters = Decenters(:,deco');

Decenters = [M1xy [1:NumSegs]' Decenters];

%% Prepare data and output
Alldata = [RES.ActNodePosListTMTM1 RES.ActNodeNumList Actstrokes];

csvwrite(actfilename,Alldata);
csvwrite(pointfilename,Pointings);
csvwrite(decfilename,Decenters);

% End of <nrSciActCSVFile>



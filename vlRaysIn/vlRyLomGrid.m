%
% File: <vlRyLomGrid.m>
%
% Syntax: vlRyLomGrid
%
% Description:
%       Initializes the list of rays to trace through the system.  The
%       format of the RAYSIN data structure is:
%           Status,Wave,i,X,Y,Z,L,M,N,D
%
%           Column 1: 	Status, only used for OUTRAYS structure, ignored on input.
%           Column 2: 	Wave, Zemax wavelength number.
%           Column 3: 	i, Ray intensity.
%           Column 4: 	X, coordinate of ray at OC.TrStart Surface
%           Column 5: 	Y, coordinate of ray at OC.TrStart Surface
%           Column 6: 	Z, coordinate of ray at OC.TrStart Surface
%           Column 7: 	L, direction cosine of the ray in the X direction
%           Column 8: 	M, direction cosine of the ray in the Y direction
%           Column 9: 	N, direction cosine of the ray in the Z direction
%           Column 10: 	D, Length of the ray, can be non-zero on input.
%
% Input Parameters:
%       none
%
% Output Parameters:
%       Creates global data structure RAYSIN
%
% Required Global Data Structures:
%       OC, RAYSIN
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
% Copyright (c) 2002,2003 Scott Roberts, National Research Council
%       Rev 1.1, Initial Release, 30 October 2002
%       Rev 1.2, Updated data structure, 13 April, 2003
%       Rev 1.3, Fixed OC.Wave bug and array init. bug, 6 May, 2002
%
% static char rcsid[] = "$Id: vlRyLomGrid.m,v 1.6 2004/07/13 19:32:14 kerleyd Exp $";
% INDENT-OFF*
% $Log: vlRyLomGrid.m,v $
% Revision 1.6  2004/07/13 19:32:14  kerleyd
% removes referance to 'currentDir', since variable undefined
%
% Revision 1.5  2004/02/19 00:28:17  kerleyd
% added gobal IM
%
% Revision 1.4  2003/12/19 22:18:34  stretchn
% Un-hard-coded array of non-missing segments
%
% Revision 1.3  2003/12/19 17:39:26  stretchn
% Removed hard-coded 28 segments
%
% Revision 1.2  2003/12/19 02:27:53  stretchn
% Now tweaks RAYSIN so that all rays hit segments
%
% Revision 1.2  2003/09/19 17:39:59  stukasa
% Header update.
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

function vlRyLomGrid

global IM;
global RAYSIN;
global OC;

% Check to make sure the OC data structure is defined, and not empty
% Must test for one of the structure members as the global declaration 
% makes exist('OC','var') true
if isempty(OC);
    error('OC Data Structure is not defined');
end

%RAYSIN = zeros(6*(OC.ConSegs-length(OC.SegMiss)),10);

Samples = OC.LOM.N;
Xoff = OC.LOM.RayOrient(1);
Yoff = OC.LOM.RayOrient(2);
Zoff = OC.LOM.RayOrient(3);
XDirCos = OC.LOM.RayOrient(4);
YDirCos = OC.LOM.RayOrient(5);
WaveNum = OC.LOM.WaveNum;

% Set up SSP, PMSP data structures
SSP = zeros(OC.NumSeqs,6);
PMSP = zeros(OC.NumNonSeqs,8);

for const = 1:6
    for seg = 1:OC.ConSegs
        row = (const - 1)*OC.ConSegs + seg;
        PMSP(row,1) = const;
        PMSP(row,2) = seg;
    end
end

%Open Zemax DDE Link
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
%     cd(currentDir);
    error(sprintf(['   IM.ZemaxPrescription = [%s]\n    Actual prescription = [%s]\n'...
            'Zemax Prescription Mismatch'],IM.ZemaxPrescription,ZemFile));
end     

% Calculate OPL0, with mirrors in the nominal position.
vlMzMoveSurfaces(SSP,PMSP,0); % set surfaces in Zemax to nominal positions

disp('    Preparing initial RAYSIN matrix...');

%Make a vector of non-missing segments
segments = [];
point = 1;
for k = 1:length(OC.SegMiss)
    segments = [segments,point:OC.SegMiss(k)-1];
    point = OC.SegMiss(k) + 1;
end
segments = [segments,point:OC.ConSegs];


[ RAYSIN,Status,Map ] = vlLmMakeRaysHs(Samples,[1:6]',segments',Xoff,Yoff,Zoff,XDirCos,YDirCos,WaveNum);
disp('    Adjusting for missed rays, calculating OPL0...');
[ OPL0, Status, IAD, IADNum, RAYSIN ] = vlLmGetWfm(RAYSIN,Samples,Map,zeros(Samples^2,1));
%Note: raysin is now slightly tweaked to account for any rays that fell
%between mirrors.

%Close zemax link
vlMzClose;

OPL0 = -OPL0; % to maintain proper sign convention

% compute OPDzero by normalizing OPL0
OPDzero = zeros(Samples^2,1);
OPDzero(Map) = max(OPL0(Map)) - OPL0(Map);

OC.LOM.OPL0 = OPL0;
OC.LOM.Map = Map;
OC.LOM.IAD = IAD;
OC.LOM.IADNum = IADNum;
disp('    Done.');
% end of vlRy1PerSeg.m

% File: TMTRCV1.m 
%
% Syntax: TMTOSMV1
%
% Description:
%       Fills OC global data structure corresponding to the optical
%       prescription for the delta design effort Oct 2006 where TMT
%       is now a RC telscope as per TMT.SEN.06.001.DRF02 Oct 13, 2006.
%       Number of segments is now 492 
%       size ( point to point ) 1.432m 
%       Seg missing is now 1,12,13,25,36,46,88,90,91
%       Rings are 13
%       Version 1 of the perscription!
%
% Input Parameters:
%       none.
%
% Output Parameters:
%       none.
%
% Required Global Data Structures:
%       OC        - Modified by this function
%
% Required Data Files:
%       none.
%       

%
% Revision History
% created by John Pazder Oct 26, 2006
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%           Herzberg Institute of Astrophysics                  %%%%%
%%%%%%      Astronomy Technology Research Group - Victoria           %%%%%
%
% (c) <2006>				        (c) <2006>
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


function TMTRCV1

global OC

%
%  
% optical parameters
OC.ZemaxUnits.lengthPerMetre=1;      % the number of zemax units per a meter.
OC.EPD = 30.100;                     % Entrance pupil diameter defined in zemax        (Zemax units)
OC.ExitPupilDist = 46.40654;           % Distance to the exit pupil from the focal plane (Zemax units)
OC.FL = 450 ;                        % System focal length      (Zemax units)
OC.ImageScale = OC.FL/206265;        % Image scale at OC.TrStop (Zemax units/arcsecond)
OC.Wave = [1.00e-6 1.64e-6 2.22e-6]; % Wavelength Values        (always meters)
OC.NumWaves = length(OC.Wave);   % # of wavelengths defined

% primary mirror segmentation
OC.NumRings = 13;             % Number of segment rings
OC.NumNonSeqs = 3*(OC.NumRings^2+OC.NumRings); % Total number of segments
OC.ConSegs = OC.NumNonSeqs/6; % Number of segments in a constellation
OC.NonSeqSegGap = 0.51/1000 ;      % Gap between primary mirror segments (Zemax units)
OC.SegMiss = [1 12 13 25 36 46 88 90 91]; % List of missing segments in each constellation 
OC.SegSize = 1.432;           % Segment Point to Point Diameter     (Zemax units)

% primary mirror shape
OC.PrimCC =  -1.00095;    % Conic Constant of the primary mirror  (0 = sphere)
OC.PrimRadCurv = 60.000;   % Radius of curvature of primary mirror (Zemax units)

% secondary mirror shape
OC.SecRadCurv = -6.227680; % Radius of curvature of primary mirror, positive for Gregorian (Zemax units) 

% mirror cell shape
OC.MirrCellCC = 0;               % Mirror support sell conic constant      (0 = sphere)
OC.MirrCellRad = OC.PrimRadCurv; % Mirror support cell radius of curvature (Zemax units)
OC.PrimOffset = OC.SegSize/3;    % Distance from mirror cell to primary    (Zemax units)

% zemax perscription indexes and names
OC.NonSeqSurf = 5;               % Zemax LDE Surface number of non-sequential surface
OC.SeqSurf = [10 19 28];          % Sequential ray trace surface #'s
OC.NumSeqs = length(OC.SeqSurf); % Number of sequential surfaces modified by the structural dynamics
OC.SeqSurfZdir = [-1 1 -1];      % Direction of surface Z axis
OC.SeqSurfName = {'Secondary' 'Tertiary' 'Focus'}; % Names of sequential surfaces 
% Specify files that define coordinate transforms of sequential surfaces
% between Z and the surface coordinate systems.
OC.ZMSeqCSTName ...
  = {'vlCsRcZsecondary_3mirr' 'vlCsRcZtertiary_3mirr' 'vlCsRcZfocus_3mirr'};
OC.TrStart = 1;                  % First trace surface
OC.TrStop = 28;                  % Last trace surface

% Initialize OC.LOM Data
OC.LOM.N = 256;                      % Number of samples
OC.LOM.RayOrient = [0 0 10.000 0 0]; % Parameters to define source (Zemax units)
OC.LOM.WaveNum = 1;                  % WaveNumber at which the LOM will be calculated
OC.LOM.LinVal = [0.50555e-3 0.50555e-3 0.50555e-3 0.00051 0.00051 0.00051 ...
        0.5055e-3 0.5055e-3 0.50555e-3 0.00051 0.00051 0.00051]; % Linearizing value (Zemax units)

% Initialize OC.FFT Data 
OC.FFT.Zp = 1;                            % Zero padding
OC.FFT.Icr = 206265*OC.Wave(OC.LOM.WaveNum)/((OC.FFT.Zp+1)* ...
  (OC.EPD/OC.ZemaxUnits.lengthPerMetre)); % Grid increment in arc seconds

% Initialize the OC.ZER Data
OC.ZER.Range = [1; 32];                   % Range of Zernikes to extract
OC.ZER.NumZer = OC.ZER.Range(2)-OC.ZER.Range(1)+1;  % Total # of Zernikes

% End of vlOpV1p0.m 
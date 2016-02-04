% File: <nrSciRawSegMotions.m>
%
% Syntax: nrSciRawSegMotions(filenamespec, outfilename)
%
% Description:
%        Saves uncorrected segment motions to a .mat file.  This data is
%        used by JPL for M1CS control system analysis.
%
% Input Parameters:
%       filenamespec - [filename] specifies pattern to match all filenames.
%           For example, '*Align0_Data.mat'.
%       outfilename - [filename] - output .mat filename specification
%
% Output Parameters:
%       The output data are the uncorrected segment motions from the TMT
%       FEA model for the case where the segments are aligned with the
%       telescope zenith pointing.  The motions are reported relative to a
%       best-fit M1CS coordinate system at each elevation angle.  
%
% Required Global Data Structures:
%         None
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
% $Id: nrSciRawSegMotions.m,v 1.1 2012/11/01 22:51:26 roberts Exp $
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


function nrSciRawSegMotions(filenamespec,outfilename)

D=dir(filenamespec);
numcases = length(D);


for ii=1:numcases
    fprintf('Found file: %s\n',D(ii).name);
    load(D(ii).name);
    
    result.(genvarname(['elev' num2str(RES.ElevAng)])).CST_M1CS_TMTM1Seg_Original=RES.CST_M1CS_TMTM1Seg_Original;
    result.(genvarname(['elev' num2str(RES.ElevAng)])).CST_M1CS_TMTM1Seg_Perturbed=RES.CST_M1CS_TMTM1Seg_Perturbed;
    result.(genvarname(['elev' num2str(RES.ElevAng)])).EulerXYZ_TMTM1Seg_orig_to_pert=RES.EulerXYZ_TMTM1Seg_orig_to_pert;
end

save(outfilename,'result');

% End of <nrSciRawSegMotions>



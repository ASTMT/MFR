% File: <vlOpPSFtoStrehl.m>, VLOT Toolbox
%
% Syntax: [strehl] = vlOpPSFtoStrehl(PSF)
%
% Description:
%       This routine computes the Strehl ratio of the image 
%       given the 'diffraction scaled' PSF of the wavefront.
%
%       The PSF is assumed to be scaled with respect to the 
%       Strehl ratio with a unit of 1 being the peak intensity
%       from a diffraction limited image.
%
%       The routine PSF_FFT returns the PSF with this scaling.
%
%       Thus, by the definition of Strehl ratio, the Strehl ratio
%       is max(PSF)
%
%       The user should understand that no interpolation is done on the
%       data and that the apropriate pupil sampling and zero pading must be
%       applied to the data in the computation of the PSF for a reasonable
%       Strehl Ratio.
%
%       This routine works with column vector formats of PSF
%
% Input Parameters:
%       PSF       - (Nx1 vector) the diffraction scaled PSF ( Point Spread Function ) in column
%                   or mxm array format.
%
% Output Parameters:
%       strehl    - (scalar) The Strehl ratio of the given PSF.
%
% Required Global Data Structures:
%       none
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
%   Copyright (c) 2003, John Pazder, National Research Council
%       Rev 1, Initial Release, June 11,2003  
%
% static char rcsid[] = "$Id: vlOpPSFtoStrehl.m,v 1.2 2003/09/19 21:22:33 stukasa Exp $";
% INDENT-OFF*
% $Log: vlOpPSFtoStrehl.m,v $
% Revision 1.2  2003/09/19 21:22:33  stukasa
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

function  strehl = vlOpPSFtoStrehl(PSF)

if  nargin ~= 1,  error('Error - Incorrect parameters, must enter one parameter, PSF') ; end

% compute Strehl ratio
strehl=max(PSF(:));

% end of function
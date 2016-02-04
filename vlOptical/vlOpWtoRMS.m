% File: <vlOpWtoRMS.m> , VLOT Toolbox
%
% Syntax: [rms] = vlOpWtoRMS(opd)
%         [rms] = vlOpWtoRMS(opd,iad)
%
%  opd      = (Nx1 vector) Optical Path difference of WF 
%  iad      = (Nx1 vector) Intensity Apodization data
%
%  pv       = (scalar) the peak valley of the wavefront
%
% Description:
%       This routine computes the RMS( Root Mean Square ) of the inputed opd.
%
%       The input data structure is either the short form, which uses the 
%       global OC data structure or the long form which requires all data 
%       to be inputed.
%
%       The units of the RMS are the same as the units of the 
%       opd inputed.
%
%       In the long form, when opd and iad are given, the opd is masked with 
%       the given iad, so the pv is only calculated for only valid aperture points. 
%
%       In the short form, when only opd is give, the opd is masked using the 
%       OC.LOM.Map, so the pv is calculated only on the valid aperture
%       points as defined by the map vector
%
%       The RMS is not flux weighted with the iad ( only masked )!!!!
%
%       Mathmaticaly, the RMS is the second moment of the sample about
%       the mean, thus RMS is given as follows:

%       RMS = sum((OPD-mean(OPD))^2)/ (number of points)
%             for all masked points on the OPD.
%
%       This routine works with column vector formats of opd and iad.
%
% Input Parameters:
%       opd       - (Nx1 vector) Optical Path difference of the wavefront at 
%                   the exit/entrance pupil. 
%                   The opd must be a column vector.
%                   The length of the iad and opd must be the same.
%
%       iad       - (Nx1 vector) The intensity apodization data of the wavefront 
%                   at the exit/entrance pupil. The data is the intensity of the
%                   wavefront at the pupil. This effectivly defines the 
%                   aperture of the telescope, but also may define atmospheric 
%                   scintillation effects. Areas of the wavefront which the
%                   wavefront is not defined have an IAD of zero.
%                   The iad must be a column vector.
%                   The length of the iad and opd must be the same.
%                   When no iad is specifed, the OC.LOM.Map vector is used to
%                   mask the opd for valid aperture points.
%
% Output Parameters:
%       rms       - (scalar) The rms ( root mean square ) of the Optical Path difference of
%                   the wavefront masked over the iad. Units are same as
%                   the units of the OPD inputed.
%
% Required Global Data Structures:
%       with one argument, OC data structure.
%       data items: OC.LOM.Map
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
% static char rcsid[] = "$Id: vlOpWtoRMS.m,v 1.6 2005/02/03 00:00:21 stretchn Exp $";
% INDENT-OFF*
% $Log: vlOpWtoRMS.m,v $
% Revision 1.6  2005/02/03 00:00:21  stretchn
% Removed unnecessary IM global struct.
%
% Revision 1.5  2004/11/08 21:28:04  dunn
% Fixed the isfield call.
%
% Revision 1.4  2003/12/19 01:13:05  stretchn
% Now masks the OPD for Zemax cases where a square OPD has been generated - these cases will be identical to the LOM at this point
%
% Revision 1.3  2003/11/24 19:46:20  stretchn
% Different RMS calculation for Zemax and LOM now.
%
% Revision 1.2  2003/09/19 22:48:56  stukasa
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


function  [rms] = vlOpWtoRMS(opd,iad)
global OC

if isfield(OC,'LOM') && isstruct(OC.LOM) && isfield(OC.LOM,'Map') && isfield(OC.LOM, 'N') && length(opd) == OC.LOM.N^2
    if nargin == 1
        % short form, use global data structure
        % mask opd with iad, extracting all values which the iad is not zero
        % using map
        opd_masked=opd(OC.LOM.Map);
    else
        % this is is non global data structure mode,
        % check for 2 input args
        if  nargin ~= 2
            error('Error - Incorrect parameters, must enter one ( opd ) or two parameters ( iad and opd )');
        end
        % check OPD and IAD are same length
        if  any( size(opd) ~= size(iad) ), error(' Error - OPD and IAD must have the same number of elements') ; end
        % mask opd with iad, extracting all values which the iad is not zero
        % use  loginal array of iad values ~=0
        opd_masked=opd(iad~=0);
    end
    
    % compute RMS via matlab standard deviation equation with flag=1 for
    % 'second moment' form of the equation
    rms=std(opd_masked,1);
else
    %All rays present will be valid if running with Zemax.  Otherwise, opd
    %will be zeros.
    rms = std(opd,1);
end
% end of function

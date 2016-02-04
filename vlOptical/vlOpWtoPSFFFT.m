% File: <vlOpWtoPSFFFT.m>, VLOT Toolbox
%
% Syntax: [PSF,icr] = vlOpWtoPSFFFT(opd)
%         [PSF,icr] = PSF_FFT(OPD,IAD,wl,zp,epd)
%         
%  opd = (Nx1 vector) optical path difference
%  iad = (Nx1 vector) intensity apodization data 
%  wl  = (scalar) wave length in microns
%  zp  = (scalar) zero padding factor
%  epd = (scalar) entrace pupil diameter in mm
%
%  PSF = (Nx1 vector) intensity distribution at image surface
%  icr = (scalar) grid increment of image surface in arcseconds.
%
% Description:
%       This routine computes the Point Spread function for the telescope
%       given the opd and iad at the entrance/exit pupil of the telescope.
%       This is the predicted intensity distribution at the image plane for
%       a point object.
%
%       The input data structure is either the short form, which uses the 
%       global OC data structure or the long form which requires all data 
%       to be inputed.
%
%       In the short form : iad=OC.LOM.IAD; zp=OC.FFT.Zp;
%       wl=OC.Wave(OC.LOM.WaveNum); epd=OC.EPD.  
%
%       The FFT ( Fast Fourier Transform ) is used to calculate the PSF.
%       The use of the FFT has some limitations and approximations. These
%       are: - the aperture must be >> then the wavelength
%            - The exit pupil distance must be >> then the wavelength.  
%            - the extended sine condition must be met. In other words, the
%              aberrations must be small enough that the mapping of the
%              entrance and exit pupils must be free of distortion.
%            - The exit pupil distance must be >> then the wavelength.  
%            - polarization is not considered
%            - a scalar approximation is sufficient, thus the F/# must be
%              slower than  1.5   
%            - The phase difference between any two points in the wavefront 
%              data is < 1/4 of a wave, or the data will phase wrap.
%
%        As a further approximation, for this routine it is assumed that the
%        point of calculation is either on axis or very close to on axis,
%        thus the actual F/ratio of the beam at the point of calculation
%        can be considered to be the F/ratio of the telescope, with the
%        cosine factor of the off axis focal length and image plane
%        curvature can be ignored. The Obliquity of projection, or tilt of
%        the image plane and the exit pupil are ignored. 
%
%        The size of the wavefront ( OPD/IAD) is assumed to be epd x epd in
%        width at the entrance pupil. Thus the opd data is equal to the epd. 
%        The scale of the psf is given by icr=206265*wl/((zp+1)*(epd)). 
%        For a circular aperture with a diameter given by epd the diameter 
%        of the diffraction limited image up to the first dark diffraction 
%        ring is given by 2.44*wl/epd. Thus the number of sample points 
%        across the PSF is given by 2.44*(zp+1). 
%
%        The total size of the PSF is given by  
%        icr*size(PSF)=206265*(size(OPD)-1)*wl/(epd). Thus for epd=21220mm, 
%        and wl=1 micron, the total size of the PSF is 2.49 arcseconds.  
%
%        For a given wl and epd, the total size of the PSF can be adjusted by
%        adjusting the number of sample points across the entrace pupil, or
%        size(OPD). The number of points should be selected such that the
%        OPD between any two points is < 0.25*wl. Increasing sampling
%        across the aperture will not increase the resultion of the PSF,
%        but will increase the size of the area sampled.  Thus for worse
%        seeing or optical quality, the sampling across the aperture should
%        be increased. If the sampling across the aperture is insufficent,
%        you may see aliasing in the PSF, image artifacts that appear to
%        wrap from one size to the other of the image.
%
%        The sampling of the PSF is increased via the zp factor. The
%        resolution across the difraction limited image is given by
%        2.44*(zp+1). The default value for zp is 1, thus 4.88 samples
%        across the diffraction limited PSF. Thus large zp values are only
%        useful for systems that are near-diffraction limited. 
%
% Input Parameters:
%       opd       - (Nx1 vector) Optical Path difference of the wavefront at the
%                   exit/entrance pupil.
%                   The opd must be a column vector.
%                   The length of the iad and opd must be the same.
%
%       iad       - (Nx1 vector) The intensity apodization data of the wavefront at the 
%                   exit/entrance pupil. The data is the intensity of the
%                   wavefront at the pupil. This effectivly defines the 
%                   aperture of the telescope, but also may define atmospheric 
%                   scintillation effects. 
%                   The iad must be a column vector.
%                   The length of the iad and opd must be the same.
%                   The default value for iad is OC.LOM.IAD 
%
%       wl        - (Scalar) Wavelength of the wavefront for the PSF calculation in microns. 
%                   The routine currently only calculates monochromatic PSF's. 
%                   The default value for the wl is OC.Wave(OC.LOM.WaveNum)
%
%       zp       -  (Scalar) Zero padding factor. This is the amount of zero padding to
%                   apply to the wavefront data before computing the PSF. The
%                   transformed wavefront ( and the PSF)  will be
%                   (zp+1)*size(OPD). Thus, zp=0, will apply no zero padding .
%                   zp must be >0 and an integer.
%                   The default value for zp is OC.FFT.Zp
%
%       epd      -  (Scalar) ( entrance pupil diameter ) The width of the opd data at the 
%                   entrance pupil of the telescope. This is not necessary the diameter
%                   of the primary mirror, as the epd may be slightly larger than the mirror.
%                   The default value for epd is OC.EPD
%
% Output Parameters:
%       PSF       - (Nx1 vector) This a column vector of the intensity distribution  
%                   at the image plane, or the point spread function for the inputed 
%                   wavefront. The size of the matrix is (zp+1)*size(OPD). The PSF is
%                   scaled such that an unaberrated image would have a peak intensity
%                   of 1. Thus, the strehl ratio is max(max(PSF)).
%
%       icr       - (Scalar) This is the grid increment at the image surface in arcseconds. 
%                   The total size of the PSF is given by icr*(size(PSF)-1). 
%                   The grid increment is given by icr=206265*wl/((zp+1)*(epd))
%                   In the short form, where only opd is specified, icr is not
%                   returned.
%
% Required Global Data Structures:
%       with one argument, OC data structure.
%       data items: OC.FFT.Zp
%                   OC.LOM.IAD
%                   OC.LOM.WaveNum
%                   OC.Wave
%                   OC.EPD
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
%       Rev 1, Initial Release, April 4,2003  
%
% static char rcsid[] = "$Id: vlOpWtoPSFFFT.m,v 1.5 2004/08/23 18:37:12 kerleyd Exp $";
% INDENT-OFF*
% $Log: vlOpWtoPSFFFT.m,v $
% Revision 1.5  2004/08/23 18:37:12  kerleyd
% addad global IM
%
% Revision 1.4  2004/08/23 16:52:12  kerleyd
% uses IM.UnitsPerMeter to do all calculations in meters
%
% Revision 1.3  2004/07/28 23:03:15  msmith
% Added icr as output parameter on simple (single-input) call.
%
% Revision 1.2  2003/09/19 22:37:22  stukasa
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

function [PSF,icr] = vlOpWtoPSFFFT(opd,iad,wl,zp,epd)

global IM

if nargin == 1
    % short form, use global data structure
    global OC
    % set global varables
    iad=OC.LOM.IAD;
    wl=OC.Wave(OC.LOM.WaveNum);
    zp=OC.FFT.Zp;
    epd=OC.EPD;
else
    % this is is non global data structure mode,
    % check for 5 input args
    if  nargin ~= 5
        error('Error - Incorrect parameters, must enter one (opd) or five parameters (opd,iad,wl,zp,epd)');
    end    
end

% check OPD and IAD are column vectors of same size
if  any( size(opd) ~= size(iad) ), error(' Error - OPD and IAD must have the same number of elements') ; end
 
% check that zp is >=0 and integer
if         zp <  0, error(' Error - zp must be > 0') ; end
if  mod(zp,1) ~= 0, error(' Error - zp must be integer') ; end 

% reshape opd and iad to square mxm
length=size(opd,1);
dim=sqrt(length);
if  mod(dim,1) ~= 0 ,error(' Error - opd and iad not properly formated, sqrt(size(opd,1)) is not an integer') ; end 
opd=reshape(opd,dim,dim);
iad=reshape(iad,dim,dim);

% compute new dim with zero padding
dim=(zp+1)*dim;

% compute the complex WF from the OPD and IAD data.
% cwf=iad*exp(i(2pi/wl)wfm)
cwf=iad.*exp((i*2*pi/wl)*(opd/IM.UnitsPerMeter));

% note naming convention, lower case is pre-transform, upper is transform
% CWF =FFT(cwf)

% compute complex amplitude distribution on image plane via FFT of CWF
% apply fftshift because of matlabs zeropoint location
% dim^2/sum(sum(IAD))is normalizing factor for ifft2 scaling.
CWF=dim^2/sum(sum(iad))*fftshift(ifft2(cwf,dim,dim)); 

% We now have the complex wavefront on image plane (amplitude and phase)
% We would like the intensity on the image plane
% this is computed via PSF=CWF*CWF'
% and reshape to column vector
PSF=reshape(CWF.*conj(CWF),dim*dim,1); 


% compute scale in arcseconds if long form
% remember wl in meters, epd is converted to meters

icr=206265*wl/((zp+1)*(epd/IM.UnitsPerMeter));

% end of function

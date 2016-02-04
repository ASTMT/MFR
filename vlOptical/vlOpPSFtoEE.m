% File: vlOpPSFtoEE.m
%
% Syntax: diameter = vlOpPSFtoEE(psf, eeFract, gridSpacing)
%

% Description:
%       Computes the diameter required to encircle the given fraction of
%       total energy within the PSF.
%
% Input Parameters:
%       psf         - (N^2 * 1) Column vector representing the PSF.
%                   
%       eeFract     - Vector containing the desired encircled energy ratios.
%                     Each element must be in the range 0 < eeFract(i) < 1.
%
%       gridSpacing - Spacing (in arcsec) between adjacent grid elements 
%                     in the PSF.  This is the "icr" output parameter from
%                     the call to vlOpWtoPSFFFT.
%
% Output Parameters:
%       radius   - Vector of diameter required to encircle the specified 
%                  energy fractions.
%
% Required Global Data Structures:
%       GLOBAL1   - (Optional Comments, listing of fields)
%       GLOBAL2   - (Optional Comments, listing of fields)
%
% Required Data Files:
%       File1.xxx - Description
%       File2.xxx - Description
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% static char rcsid[] = "$Id: vlOpPSFtoEE.m,v 1.5 2004/04/02 23:56:59 msmith Exp $";
% INDENT-OFF*
% $Log: vlOpPSFtoEE.m,v $
% Revision 1.5  2004/04/02 23:56:59  msmith
% Added line to initialize lo in case loop condition fails on first test
% (ie. loop does not execute).
%
% Revision 1.4  2004/03/11 20:15:07  msmith
% Added gridSpacing as a parameter and applied scaleFactor.
% Converted return values from radii to diameters.
%
% Revision 1.3  2004/03/10 04:24:57  msmith
% Converted eeFract from scalar to vector.
%
% Revision 1.2  2004/03/08 22:50:00  msmith
% Removed x and y as input parameters (always computes centroid now).
% Store distances (not squares) in dist array.
%
% Revision 1.1  2004/02/02 04:18:52  msmith
% Initial revision
%
% INDENT-ON*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%           Herzberg Institute of Astrophysics                  %%%%%
%%%%%%      Astronomy Technology Research Group - Victoria           %%%%%
%
% (c) 2004				            (c) 2004
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

function diameter = vlOpPSFtoEE(psf, eeFract, gridSpacing)

if nargin < 1
    error('Error - vlOpPSFtoEE(psf, eeFract) requires at least one parameter');
end

if nargin < 2
    eeFract = 0.5;
end

if min(eeFract) <= 0.0
    error('Error - eeFract values must be greater than 0.0');
end

if max(eeFract) >= 1.0
    error('Error - eeFract values must be less than 1.0');
end

global OC;

% Create two-dimensional arrays

len = prod(size(psf));
dim = sqrt(len);

dist = zeros(dim, dim);


% Compute X and Y centroids  (which index is X?  row or column?)
% This should probably be done using vector operations...

psf = reshape(psf, dim, dim);

sx = 0.0;
sy = 0.0;
ss = 0.0;

for col=1:dim
  for row=1:dim
    sx = sx + psf(row,col) * col;
    sy = sy + psf(row,col) * row;
    ss = ss + psf(row,col);
  end
end


% Divide the centroid sums by the number of grid points

x = sx / ss;
y = sy / ss;


% Compute distance of each grid point from target position (x,y)

for col=1:dim
    dx = col - x;
    dx2 = dx * dx;
    for row=1:dim
        dy = row - y;
        dy2 = dy * dy;
        dist(row,col) = sqrt(dx2 + dy2);
    end
end


% Reshape the distance values as a column vector so that it can be joined
% to the psf column vector.

psf = reshape(psf, dim*dim, 1);

dist = reshape(dist, dim*dim, 1);

% Create a dim^2 * 2 array containing the sorted distances in the first
% column and the corresponding energy values in the second column.

tmp = sortrows([ dist, psf ]);

diameter = zeros(length(eeFract), 1);

sumPSF = sum(psf);

% Wrap the following conversions into one scale factor:
% - Convert from radius to diameter (2)
% - Convert from grid spacing to arcsec
% - Convert from arcsec to mm in focal plane

scaleFactor = 2.0 * gridSpacing * OC.ImageScale;

for i=1:length(eeFract)

    % Determine the target amount of energy to include in the circle.

    targetEE = eeFract(i) * sumPSF;


    % Initialize the cumulative encircled energy sum to the closest grid point.

    ee = tmp(1,2);
    eePrev = 0.0;
    ndx = 1;
    lo = tmp(ndx, 1);


    % Add one grid element at a time until the encircled energy meets or
    % exceeds the target value.

    while (ee < targetEE) && (ndx < len)
        eePrev = ee;
        lo = tmp(ndx, 1);
        ndx = ndx + 1;
        ee = ee + tmp(ndx, 2);
    %    if (ndx < 11)
    %    disp(sprintf('%d  %10.4f  %10.4f  %10.6f / %10.6f', ndx, tmp(ndx,1), tmp(ndx,2), ee, targetEE));
    %    end;
    end

    % Displaying or returning these values can be useful for debugging.

    hi = tmp(ndx, 1);

    % Linearly interpolate the radius value that would have the target EE.
    radius = lo + (targetEE - eePrev) * (hi - lo) / (ee - eePrev);

    diameter(i) = radius * scaleFactor;

end  % for i

% end of function

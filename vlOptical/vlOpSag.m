% File: <vlOpSag.m>, NRCIM Toolbox
%
% Syntax: [z dzdy] = vlOpSag(c,bs,y)
%
% Description:
%       This routine calculates the sag of a conic surface.
%   
%       Reference: "Reflecting Telescope Optics I", R.N. Wilson, page 80,
%       equation 3.75.  Also see the Standard surface type in the Zemax manual.
%       The derivative was calculated symbolically in MathCAD.
%
%                   z = c*y^2/(1+sqrt(1-(1+bs)*c^2*y^2))
%
%                   dz/dy = c*y/sqrt(1-(1+bs)*c^2*y^2)
%
% Input Parameters:
%       c         - curvature of the surface at the vertex (1/radius)
%       bs        - conic constant (parabola = -1)
%                      bs < -1 for hyperbolas
%                      bs = -1 for parabolas
%                      -1 < bs < 0 for ellipses
%                      bs = 0 for spheres
%                      bs > 0 for oblate ellipses%
%       y         - N*1 column vector of distances from the vertex in the 
%                   plane of the tangent to the vertex
%
% Output Parameters:
%       z    - N*1 column vector of departures from the plane of the 
%              tangent to the surface vertex
%      dydz  - N*1 column vector of derivatives
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
%   Copyright (c) 2002, Scott Roberts, National Research Council
%       Rev 1.1, Initial Release, 24 October 2002  
%
% static char rcsid[] = "$Id: vlOpSag.m,v 1.4 2012/10/24 15:56:09 roberts Exp $";
% INDENT-OFF*
% $Log: vlOpSag.m,v $
% Revision 1.4  2012/10/24 15:56:09  roberts
% Updated with simplified equations
%
% Revision 1.3  2006/06/05 16:47:10  msmith
% Converted to vectorized form to allow y parameter to be a vector.
%
% Revision 1.2  2003/09/19 21:50:53  stukasa
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

function [z, dzdy] = vlOpSag(c, bs, y)

A = sqrt(1-(1+bs)*c^2*y.^2);

z = c*y.^2 ./ (1+A);

dzdy = c*y ./ A;

% End of vlOpSag.m

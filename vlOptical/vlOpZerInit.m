% File: vlOpZerInit.m
%
% Syntax: vlOpZerInit
%
% Description:
%
%   This routine initialize all the values needed for the computation of
%   the coefficients of the Zernike polynomials.
%
%   First, it calculates the average outer radius of the pupil. This value
%   is needed for the calculation of the zernike polynomials.
%
%   Then, it calculates the normalized Zernike polynomials and stores it in
%   OC.ZER.Poly. Finally, it calculates the inverse gmm matrix witch is
%   needed to do the appropriate correction to the coefficients in order to
%   have a linearly independent set of coefficients. It stores this matrix
%   in OC.ZER.InvGmm.
%
% Input Parameters:
%
% Output Parameters:
%
% Required Global Data Structures:
%       OC.LOM.N
%       OC.LOM.IAD
%       OC.LOM.IADNum
%       OC.ZER.Range
%
% Required Files:
%   vlOpGetZer.m
%   vlOpZerNumero.m
%   vlOpScalarProd.m
%   vlOpGmm.m
%

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% static char rcsid[] = "$Id: vlOpZerInit.m,v 1.9 2008/07/31 15:41:07 msmith Exp $";
% INDENT-OFF*
% $Log: vlOpZerInit.m,v $
% Revision 1.9  2008/07/31 15:41:07  msmith
% Set start of Zernike range to 1 in call to vlOpGetZer.
%
% Revision 1.8  2005/02/01 19:36:58  msmith
% Removed getRadius function.  Set diameter to OC.LOM.N-1.
%
% Revision 1.7  2005/01/11 20:38:38  msmith
% Moved radius estimation code into separate function.
% Set diameter to OC.LOM.N-1 to make consistent with LOM.
%
% Revision 1.6  2003/09/19 23:04:37  stukasa
% Header update.
%
% Revision 1.5  2003/08/21 19:11:00  msmith
% Changed bitwise & and | to logical && and ||.
%
% Revision 1.4  2003/08/10 15:36:39  vlotim
% *** empty log message ***
%
% Revision 1.3  2003/06/24 21:36:04  lavigne
% Header changed
%
% INDENT-ON*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%           Herzberg Institute of Astrophysics                  %%%%%
%%%%%%      Astronomy Technology Research Group - Victoria           %%%%%
%
% (c) <2003-2008>				    (c) <2003-2008>
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


function vlOpZerInit

global OC

N = OC.LOM.N;

IAD=reshape(OC.LOM.IAD,N,N);

d = N - 1;

% The program getzer_norm is called to calculate the required Zernike
% polynomials.  Always start at Zernike 1.
OC.ZER.Poly=vlOpGetZer(1,OC.ZER.Range(2),N,d,IAD,OC.LOM.IADNum);

% The Gmm matrix is computed and stored in OC.ZER.Gmmmm for later use.
OC.ZER.Gmm=vlOpGmm(OC.ZER.Poly,IAD,OC.LOM.IADNum,OC.ZER.Range(2));

% The inverse of Gmm is calculated
OC.ZER.InvGmm=inv(OC.ZER.Gmm);

% End of vlOpZerInit.m

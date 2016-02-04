% File: vlOpWtoZer.m
%
% Syntax: c=vlOpWtoZer(phi)
%
% Description:
%           This function transforms the OPD vector into an OPD matrix and it
%       projects it on each defined Zernike polynomial. A first set of
%       coefficients are then found. A correction is needed since the
%       polynomials are not necessarily orthogonal on the aperture. This
%       correction is made by multiplying those coefficients by the inverse
%       Gmm matrix containing the dependance of each polynomials in terms of
%       the others.
%
% Input Parameters:
%       phi       - (Oc.LOM.N^2 x 1) OPD vector for which the Zernike polynomials 
%                   coefficients are calculated
%
% Output Parameters:
%       c         - (OC.ZER.Range(2) x 1) linearly independent Zernike coefficients
%
% Required Global Data Structures:
%       OC.LOM.N
%       OC.ZER.Range
%       OC.ZER.Poly
%       OC.ZER.IAD
%       OC.ZER.IADNum
%
% Required Files:
%       vlOpScalarProd.m
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% static char rcsid[] = "$Id: vlOpWtoZer.m,v 1.6 2004/04/19 23:58:14 msmith Exp $";
% INDENT-OFF*
% $Log: vlOpWtoZer.m,v $
% Revision 1.6  2004/04/19 23:58:14  msmith
% Convert Zernike lookup value to doubles when passing to vlOpScalarProd.
%
% Revision 1.5  2003/12/19 01:12:33  stretchn
% Now allows the calculation for Zemax
%
% Revision 1.4  2003/12/05 19:57:23  stretchn
% Added check that OC.OpticalEngine = LOM
%
% Revision 1.3  2003/11/25 21:39:10  stretchn
% Only calculate if OE is LOM
%
% Revision 1.2  2003/06/26 19:24:33  lavigne
% Permission changed
%
% Revision 1.1  2003/06/25 19:50:12  lavigne
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

function c=vlOpWtoZer(phi)

global OC

IAD=reshape(OC.LOM.IAD,OC.LOM.N,OC.LOM.N);

% transform the OPD vector into an OPD matrix
phi2=reshape(phi,OC.LOM.N,OC.LOM.N);

% Project the OPD matrix on each circular polynomial to find the
% coefficients
w=1:OC.ZER.Range(2);
a=zeros(1,OC.ZER.Range(2));

for i=1:OC.ZER.Range(2)
    a(i)=vlOpScalarProd(phi2, double(OC.ZER.Poly(:,:,w(i))), IAD, OC.LOM.IADNum);
end

% Correct the coefficients to have a set of linearly independant
% coefficients
c=OC.ZER.InvGmm*a';

% End of vlOpWtoZer.m
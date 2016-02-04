% File: vlOpGmm.m
%
% Syntax: gmm=vlOpGmm(zern,pup,tpup,n)
%
% Description:
%   This function calculates the matrix containing the scalar
%   product of the i_th zernike polynomial by the j_th zernike polynomial
%   The inverse of the resultant matrix allows the user to make a correction to the
%   calculated coefficients so that they become linearly independant.
%
% Input Parameters:
%    zern - (OC.LOM.N x OC.LOM.N x OC.ZER.Range(2)) Represents the
%           value of a given Zernike's polynomial for each ray composing 
%           the wavefront deformation matrix (OC.ZER.Poly).
%     pup - (OC.LOM.N^2 x 1) The pupil matrix (OC.LOM.IAD)
%    tpup - (scalar) Number of valid rays on the pupil (OC.LOM.IADNum)
%       n - (scalar) Number of the maximal Zernike polynomial to calculate
%
% Output Parameters:
%     Gmm - matrix containing the scalar product of each
%           polynomial by each other polynomial (OC.ZER.Gmm)
%
% Required Global Data Structures:
%
% Required files:
%   vlOpScalarProd.m
%

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% static char rcsid[] = "$Id: vlOpGmm.m,v 1.3 2004/04/19 23:57:07 msmith Exp $";
% INDENT-OFF*
% $Log: vlOpGmm.m,v $
% Revision 1.3  2004/04/19 23:57:07  msmith
% Changes required by converting Zernike lookup table to singles.
%
% Revision 1.2  2003/06/26 19:20:59  lavigne
% Permission changed
%
% Revision 1.1  2003/06/25 19:48:34  lavigne
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


function gmm=vlOpGmm(zern,pup,tpup,n);

% initialize a position vector
w=1:n;

% initialize the gmm matrix
gmm=zeros(n);

% calculate the cross product of every polynomials
for i=1:n
    zi = double(zern(:,:,w(i)));
    for j=1:i
        zj = double(zern(:,:,w(j)));
        gmm(i,j)=vlOpScalarProd(zi, zj, pup, tpup);
        if i~=j
            gmm(j,i)=gmm(i,j);
        end
    end
end

% End of vlOpGmm.m
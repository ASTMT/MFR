% File: vlOpScalarProd.m
%
% Syntax: c=vlOpScalarProd(mode1,mode2,pup,tpup)
%
% Description:
%       This routine computes the scalar product of one matrix by another
%       of the same size over the pupil defined.
%
% Input Parameters:
%       mode 1   - (OC.LOM.N x OC.LOM.N) First matrix to be multiplied (OC.ZER.Poly(:,:,i))
%       mode 2   - (OC.LOM.N x OC.LOM.N) Second matrix to be multiplied (OC.ZER.Poly(:,:,j))
%          pup   - (OC.LOM.N x OC.LOM.N) Pupil defined (OC.LOM.IAD)
%         tpup   - (scalar) Number of valid pixels in the pupil (OC.LOM.IADNum)
%
% Output Parameters:
%          scp   - (scalar) Value of the scalar product
%
% Required Global Data Structures:
%
% Required Data Files:
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% static char rcsid[] = "$Id: vlOpScalarProd.m,v 1.2 2003/06/24 21:34:14 lavigne Exp $";
% INDENT-OFF*
% $Log: vlOpScalarProd.m,v $
% Revision 1.2  2003/06/24 21:34:14  lavigne
% Header changed
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

function scp=vlOpScalarProd(mode1,mode2,pup,tpup)

scp=sum(sum(mode1.*mode2.*pup))/tpup;

% End of vlOpScalarProd.m
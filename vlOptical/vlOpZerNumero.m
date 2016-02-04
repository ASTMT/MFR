% File: vlOpZerNumero.m
%
% Syntax: num=vlOpZerNumero(zn)
%
% Description:
%       This routine calculates the value of the integers n and m for the  
%       zn_ith Zernike polynomial.
%
% Input Parameters:
%           zn    - (scalar) number of the Zernike polynomial for which n and m has to be
%                   calculated
%
% Output Parameters:
%          num    - (1x2 array) Contains the value of n in num(1) and the value
%                   of m in num(2).
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
% static char rcsid[] = "$Id: vlOpZerNumero.m,v 1.3 2003/08/10 23:41:00 msmith Exp $";
% INDENT-OFF*
% $Log: vlOpZerNumero.m,v $
% Revision 1.3  2003/08/10 23:41:00  msmith
% Changed the nested loops for computing n and m into closed form
% expressions.
%
% Revision 1.2  2003/06/24 21:33:54  lavigne
% Headers changed
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


function num=vlOpZerNumero(zn)

nq = floor(0.5 * sqrt(8 * zn - 7) - 0.5);
mq = zn - 1 - 0.5 * nq * (nq+1);
if mod(nq,2) ~= mod(mq,2)
  mq = mq + 1;
end

num = [ nq, mq ];
% End of vlOpZerNumero.m

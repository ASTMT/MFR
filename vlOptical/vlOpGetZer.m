% File: vlOpGetZer.m
%
% Syntax: zern=vlOpGetZer(first,zn,dim,d,pup,tpup)
%
% Description:
%     This routine computes a 3D matrix in which are stored the normalized zernike's
%   polynomials. 
%   The two first dimensions give the value of a given
%   zernike's polynomial for each rays composing the wavefront perturbation
%   matrix.
%     This routine should not take more than a few minutes to operate even for
%   a large number of polynomials.
%     The algorithm is largely inspired by an IDL function of the same name
%   written by Jean-Pierre Veran.
%
% Input Parameters:
%   first - (scalar) OC.LOM.Range(1)
%    last - (scalar) OC.LON.Range(2)
%     dim - (scalar) OC.LOM.N
%       d - (scalar) The pupil diameter
%     pup - (Oc.LOM.N^2 x 1) OC.LOM.IAD
%    tpup - (scalar) OC.LOM.IADNum
%
%
% Output Parameters:
%    zern - (Oc.LOM.N x OC.LOM.N x OC.ZER.Range(2)) OC.ZER.Poly
%
% Required Global Data Structures:
%
% Required files:
%   vlOpZerNumero.m
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% static char rcsid[] = "$Id: vlOpGetZer.m,v 1.6 2005/02/01 19:35:27 msmith Exp $";
% INDENT-OFF*
% $Log: vlOpGetZer.m,v $
% Revision 1.6  2005/02/01 19:35:27  msmith
% Changed x and y values from i-R to R-i.
%
% Revision 1.5  2004/04/19 23:57:53  msmith
% Converted Zernike lookup table to singles.
%
% Revision 1.4  2003/07/07 17:24:57  lavigne
% Correction to the normalization
%
% Revision 1.3  2003/06/24 21:35:02  lavigne
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


function zern=vlOpGetZer(first,zn,dim,d,pup,tpup)


% initialize the matrix in which the zernike's polynomials will be stored
zern = repmat(single(0.0), [ dim, dim, zn ]);

R=(dim-1)/2;


% find the distance from the center of the pupil and the angle between the
% radius and an horizontal line going through the center of the aperture
x=zeros(dim);
y=zeros(dim);

for i=0:dim-1
    x(:,i+1)=double(R-i);
    y(i+1, :) = double(R-i);
end


r=sqrt(x.^2+y.^2)*2/d;
theta=atan2(y,x);

% Calculate the circular Zernike polynomials
for j=first:zn
    % find the value of n and m for each polynomials
    znm=vlOpZerNumero(j);
    n=double(znm(1));
    m=double(znm(2));
    z=zeros(dim);
    
    % If the number of the polynomial is odd then m is negative and the
    % correction due to the angle has to be a sine
    if mod(j,2)==1
        for i=0:(n-m)/2
            z=double(z+(-1)^i*r.^(n-2*i)*gamma(n-i+1)/(gamma(i+1)*gamma((n+m)/2-i+1)*gamma((n-m)/2-i+1)));
        end
        if m==0
            z=z*sqrt(n+1);
        else
            z=z*sqrt(2*(n+1)).*sin(m*theta);                        
        end
        
    % If the number of the polynomial is even then m is positive and the
    % correction due to the angle has to be a cosine
    else
        for i=0:(n-m)/2
            z=double(z+(-1)^i*r.^(n-2*i)*gamma(n-i+1)/(gamma(i+1)*gamma((n+m)/2-i+1)*gamma((n-m)/2-i+1)));
        end
        if m==0
            z=z*sqrt(n+1);
        else
            z=z*sqrt(2*(n+1)).*cos(m*theta);
        end
    end
    
    % normalize the value obtained so that the rms of each polynomial has
    % the unit value
    z_temp=z.*pup;
    if n==0 & m==0
        alpha2=sum(sum(z_temp.^2))/tpup;
    else
        m=sum(sum(z_temp))/tpup;    
        alpha2=sum(sum(z_temp.^2))/tpup-m^2;
        z=z-m;
    end
    z=z/sqrt(alpha2);
    zern(:,:,j) = single(z);
end

% End of vlOpGetZer.m

% File: <vlAnFormatInt.m>
%
% Syntax: vlAnFormatInt(filename)
%
% Description:
%       Reads in a file of integer numbers, formats to 
%		set # of spaces per number, as per Ansys file read 
%		format requirements, & re-writes file.
%
% Input Parameters:
%       filename  - name of file containing integer numbers to be formatted
%
% Output Parameters:
%       writes correctly formatted numbers to a new file
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
% static char rcsid[] = "$Id: vlAnFormatInt.m,v 1.2 2003/09/25 21:38:11 stukasa Exp $";
% INDENT-OFF*
% $Log: vlAnFormatInt.m,v $
% Revision 1.2  2003/09/25 21:38:11  stukasa
% Added header.
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



function vlAnFormatInt(filename)

f = load(filename);
[rows,cols] = size(f);

fid = fopen(filename,'w');

for ii = 1 : rows
    for jj = 1 : cols
        fprintf(fid,'%10d',f(ii,jj));
    end    
    fprintf(fid,'\r\n');
end 
fclose(fid);
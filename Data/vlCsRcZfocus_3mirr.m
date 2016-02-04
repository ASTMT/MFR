% File: vlCsRcZfocus_3mirr.m
%
% Syntax: T = vlCsRcZfocus_3mirr
%
% Description:
%       Returns a 3x4 CST transform matrix, representing a co-ordinate system 
%       transform M1CS to focal plane.
%
% Input Parameters:
%       None
%
% Output Parameters:
%       T = 3x4 CST matrix representing transformation
%
% Required Global Data Structures:
%       CST - but this does not appear to be used
%
% Required Data Files:
%       None
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: vlCsRcZfocus_3mirr.m,v 1.1 2006/11/18 00:54:08 msmith Exp $
%
% INDENT-OFF*
% $Log: vlCsRcZfocus_3mirr.m,v $
% Revision 1.1  2006/11/18 00:54:08  msmith
% Initial revision
%
%
% INDENT-ON*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%           Herzberg Institute of Astrophysics                  %%%%%
%%%%%%      Astronomy Technology Research Group - Victoria           %%%%%
%
% (c) <2006>                              (c) <2006>
% National Research Council               Conseil national de recherches
% Ottawa, Canada, K1A 0R6                 Ottawa, Canada, K1A 0R6
% All rights reserved                     Tous droits reserves
%                                    
% NRC disclaims any warranties,           Le CNRC denie toute garantie
% expressed, implied, or statu-           enoncee, implicite ou legale,
% tory, of any kind with respect          de quelque nature que se soit,
% to the software, including              concernant le logiciel, y com-
% without limitation any war-             pris sans restriction toute
% ranty of merchantability or             garantie de valeur marchande
% fitness for a particular pur-           ou de pertinence pour un usage
% pose.  NRC shall not be liable          particulier.  Le CNRC ne
% in any event for any damages,           pourra en aucun cas etre tenu
% whether direct or indirect,             responsable de tout dommage,
% special or general, consequen-          direct ou indirect, particul-
% tial or incidental, arising             ier ou general, accessoire ou
% from the use of the software.           fortuit, resultant de l'utili-
%                                         sation du logiciel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = vlCsRcZfocus_3mirr

dist = vlCsRcGetDistances_3mirr;

T = vlCsMult(vlCsTrans([0 0 dist.M1M3]'), ...
             vlCsRotY(-90), vlCsRotZ(90), vlCsTrans([0 0 -dist.M3Focus]'));

% End of vlCsRcZfocus_3mirr

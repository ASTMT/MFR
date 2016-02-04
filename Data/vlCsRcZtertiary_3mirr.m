% File: vlCsRcZtertiary_3mirr.m
%
% Syntax: T = vlCsRcZtertiary_3mirr
%
% Description:
%       Returns a 3x4 CST transform matrix, which transforms a point from 
%       the M3 coordinate system to the primary mirror coordinate system.
%
%       T translates P_M3 from the M3 CS to the M1 CS of the IM.
%       P_M1 = vlCsPMult(T, P_M3)
%
% Input Parameters:
%      None
%
% Output Parameters:
%       T = 3x4 CST matrix representing transformation
%
% Required Global Data Structures:
%       CST - but this does not appear to be used
%
% Required Data Files:
%       vlCsRcGetDistances_3mirr.m
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: vlCsRcZtertiary_3mirr.m,v 1.1 2006/11/18 00:53:52 msmith Exp $
%
% INDENT-OFF*
% $Log: vlCsRcZtertiary_3mirr.m,v $
% Revision 1.1  2006/11/18 00:53:52  msmith
% Initial revision
%
%
% INDENT-ON*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%           Herzberg Institute of Astrophysics                  %%%%%
%%%%%%      Astronomy Technology Research Group - Victoria           %%%%%
%
% (c) <2006>                             (c) <2006>
% National Research Council              Conseil national de recherches
% Ottawa, Canada, K1A 0R6                Ottawa, Canada, K1A 0R6
% All rights reserved                    Tous droits reserves
%                               
% NRC disclaims any warranties,          Le CNRC denie toute garantie
% expressed, implied, or statu-          enoncee, implicite ou legale,
% tory, of any kind with respect         de quelque nature que se soit,
% to the software, including             concernant le logiciel, y com-
% without limitation any war-            pris sans restriction toute
% ranty of merchantability or            garantie de valeur marchande
% fitness for a particular pur-          ou de pertinence pour un usage
% pose.  NRC shall not be liable         particulier.  Le CNRC ne
% in any event for any damages,          pourra en aucun cas etre tenu
% whether direct or indirect,            responsable de tout dommage,
% special or general, consequen-         direct ou indirect, particul-
% tial or incidental, arising            ier ou general, accessoire ou
% from the use of the software.          fortuit, resultant de l'utili-
%                                        sation du logiciel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = vlCsRcZtertiary_3mirr

dist = vlCsRcGetDistances_3mirr;

T = vlCsMult(vlCsTrans([0 0 dist.M1M3]'),vlCsRotY(180),vlCsRotZ(-90),vlCsRotX(45))

% End of vlCsRcZtertiary_3mirr

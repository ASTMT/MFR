% File: vlSciPlotActStroke.m
%
% Syntax: h = vlSciPlotActStroke(pos,act,plotTitle)
%
% Description:
%       Creates a plot of actuator stroke requirements.  Note that the
%       actuator positions required for this function can be generated with
%       vlSciGetActNodes.
%
% Input Parameters:
%       pos - array of actuator positions in the TMT M1 coordinate system
%               [length(RES.YNsActNlsSpiral) x 3]
%       act - vector of actuator positions or stroke requirements
%               [length(RES.YNsActNlsSpiral) x 3]
%       plotTitle - Plot Title [string]
%
% Output Parameters:
%       h - handle to the figure.
%
% Required Global Data Structures:
%       none
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: vlSciPlotActPos.m,v 1.1 2006/02/22 17:29:42 roberts Exp $
%
% INDENT-OFF*
% $Log: vlSciPlotActPos.m,v $
% Revision 1.1  2006/02/22 17:29:42  roberts
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
function h = vlSciPlotActStroke(pos,act,plotTitle)

figure;
%h = scatter(pos(:,1),pos(:,2),5,RES.YNsActNlsSpiral,'filled');
h = scatter(pos(:,1),pos(:,2),15,act,'filled');
colorbar;
title(plotTitle);

% End of vlSciPlotActPos
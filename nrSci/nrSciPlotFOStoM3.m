% File: <nrSciPlotFOStoM3.m>
%
% Syntax: asets = nrSciPlotFOStoM3(filenamespec)
%
% Description:
%       Plots the position of the displaced M3 mirror from the displaced
%       focal plane coordinate system (non-rotating with elevation angle)
%       over a set of elevation angles.  Also returns the angle sets.
%
% Input Parameters:
%       filenamespec - (string) specifies pattern to match all filenames.
%           For example, '*Align0_Data.mat'.
%
% Output Parameters:
%       asets - Rows of elevation angle, then displacements and EulerXYZ
%       angles for each point plotted.
%
% Required Global Data Structures:
%       None - structures are loaded from input files
%
% Required Data Files:
%       Output Data files from MFR runs of gravity loading analysis 
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: nrSciPlotFOStoM3.m,v 1.1 2012/10/31 16:42:12 roberts Exp $
%
% INDENT-OFF*
% $Log: nrSciPlotFOStoM3.m,v $
% Revision 1.1  2012/10/31 16:42:12  roberts
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


function asets = nrSciPlotFOStoM3(filenamespec)

D=dir(filenamespec);

text_x_off = -1.5e-4;
text_x_mult = -1.5e-4;
text_y_off = -1e-5;
text_y_mult = 3e-4;
fig_size_x = 800;
fig_size_y = 600;

asets = [];

for ii=1:length(D)
    fprintf('Found file: %s\n',D(ii).name);
    elevang = vlUtGetIntFromString(D(ii).name);
    load(D(ii).name);
    T_FOSA_M3A = vlCsMult(vlCsInv(GRAVRES.T_M1A_FOSA),GRAVRES.T_M1A_M3A);
    %fprintf('%d,%.3e,%.3e,%.6f,%.3e,%.3e,%.3e\n',elevang,GRAVRES.EXYZ_T_FOSDP_NoElev_M3DP);
    asets = [asets; [elevang GRAVRES.EXYZ_T_FOSDP_NoElev_M3DP]];
end

asets=sortrows(asets);
plot(asets(:,2),asets(:,3),':ob','MarkerFaceColor','b','MarkerSize',3,'LineWidth',2)
pos=get(gcf,'Position');
set(gcf,'Position',[pos(1)-(fig_size_x-pos(3)/2) pos(2)-(fig_size_y-pos(4)/2) fig_size_x fig_size_y]);
axis fill;
 
x = asets(:,2);
y = asets(:,3);
a = asets(:,1); b = num2str(a); c = cellstr(b);
sine = -1;
for ii = 1:length(D)
    dx(ii) = x(ii) + text_x_off + sine * text_x_mult;    
    dy(ii) = y(ii) + text_y_off + sine * text_y_mult;  
    sine = sine*-1;
end

text(dx, dy, c);
%str = ['Transform FOS to M3: translation then EulerXYZ angles (dx,dy,dz,rx,ry,rz) = ' num2str(GRAVRES.EXYZ_T_FOSDP_NoElev_M3DP)];
title({'Position of displaced M3 in displaced Focal Plane Coordinates' ; 'Focal Plane coordinates do not rotate with elevation angle';'Elevation angle in [deg]'});
xlabel('Position in Focal Plane CS X [m]');
ylabel('Position in Focal Plane CS Y [m]');

saveas(gcf,'FOStoM3Position');
saveas(gcf,'FOStoM3Position','tif');

fprintf('\n\nAngle Sets for transform from FOS displaced to M3 displaced\n');
fprintf('\tEl Ang,\tDisplacementXYZ,\t\t\t EulerXYZ angles\n');
disp(asets);

% End of <nrSciPlotFOStoM3>



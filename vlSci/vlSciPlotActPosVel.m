% File: vlSciPlotActPosVel.m
%
% Syntax: [h] = h=vlSciPlotActPosVel(act,pos, plottitle)
%
% Description:
% vlSciPlotActPosVel(act,pos,'title') takes as an input a grid of actuator
% deflections for respective angles as an M x N matrix - act, the input 
% angles at which these deflections took place as an N x 1 vector, and a
% plot title, for which circumstances the position analysis is outputted.
% were produced. For the maximum and minimum displacement graph, the
% actuator numbers are given as well to allow the user to know if an
% individual actuator is showing up multiple times. For both position and 
% "velocity", where velocity is the position change of the actuator per
% change in elevation angle; 25th, 50th and 75th percentile along with mean
% are plotted as subplots on a single plot. For each point which satisfies 
% a criteria, the entire history of that actuator is plotted. The
% "velocity" is calculated by quadratic fit, with three points taken at a
% time; and the value is always that of the slope at the middle of the mini
% parabola, except for the end points where the endvalues of slope are
% taken as compromise. This function requires at least 3 components in the
% position vector.
%
% Input Parameters:
%       act       -  (M x N) matrix
%                    M is number of actuators, N is the number of elevation
%                    angles
%       pos       -  (Nx1) vector
%                    position angle
%       plottitle -  this is a string and will be the title and the
%       maxima/minima graphs of the position and velocity information
%     
%
% Output Parameters:
%       h         -  plot handle
%

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: vlSciPlotActPosVel.m,v 1.2 2006/06/07 13:31:05 roberts Exp $
%
% INDENT-OFF*
% $Log: vlSciPlotActPosVel.m,v $
% Revision 1.2  2006/06/07 13:31:05  roberts
% all plots now have same scale
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





function h=vlSciPlotActPosVel(act,pos, plottitle)
%vlSciPlotActPosVel(act,pos,'title') takes as an input a grid of actuator
%deflections for respective angles as an M x N matrix - act, the input 
%angles at which these deflections took place as an N x 1 vector, and a
%plot title, for which circumstances the position analysis is outputted.
%were produced. For the maximum and minimum displacement graph, the
%actuator numbers are given as well to allow the user to know if an
%individual actuator is showing up multiple times. 

r=size(act,1);
c=size(act,2);
                
% below, the position analysis begins.
figure

maximal_pos_indices=cell(1,c);
minimal_pos_indices=cell(1,c);
maximal_pos_array=zeros(c,c);
minimal_pos_array=zeros(c,c);

%this subplot is for the minima and maxima
for i = 1:c
    [maxi, ndx] = max(act(:,i));
    maximal_pos_array(:,i)=act(ndx,:);
    maximal_pos_indices{i}=sprintf('max ar angle %.0f, %.0f',pos(i), ndx);
    [mini, ndy] = min(act(:,i));
    minimal_pos_array(:,i)=act(ndy,:);
    minimal_pos_indices{i}=sprintf('min at angle %.0f, %.0f',pos(i), ndy);
end

% the line below is meant to concatenate the matrices so
% that they fit in as a single argument in the legend.
legendarray=[maximal_pos_indices minimal_pos_indices];
plot(pos, maximal_pos_array, '-o', pos, minimal_pos_array, '-o')
legend(legendarray)

xlabel ('elevation (degrees)')
ylabel('actuator displacement (m)')
title(plottitle)

figure
hold on

sortedlist = sort(act);

%this subplot will be for the 75th percentiles
% indices75=cell(1,c);
array75=zeros(c,c);
percentile75=zeros(1,c);
for i = 1:c
    value75=sortedlist(ceil(0.75*r),i);
    k=find(act(:,i)==value75);
    array75(:,i)=act(k(1),:);
    % indices75{i}=sprintf('actuator %.0f',ndx);
    % plot(pos,act(k(1),:),'o-')
end
average = mean(act);

%this subplot will be for the 50th percentiles, it rounds upward
array50=zeros(c,c);
percentile50=zeros(1,c);
for i = 1:c
    value50=sortedlist(ceil(0.50*r),i);
    k=find(act(:,i)==value50);
    array50(:,i)=act(k(1),:);
    % plot(pos,act(k(1),:),'o-')
end

%this subplot will be for the 25th percentiles
array25=zeros(c,c);
percentile25=zeros(1,c);
for i = 1:c
    value25=sortedlist(floor(0.25*r),i);
    k=find(act(:,i)==value25);
    array25(:,i)=act(k(1),:);
    % plot(pos,act(k(1),:),'o-')
end

UpperYLimit = max([max(max(average)) max(max(array75)) max(max(array50)) max(max(array25))]);
LowerYLimit = min([min(min(average)) min(min(array75)) min(min(array50)) min(min(array25))]);

subplot(2,2,1)
plot(pos, average,'o-')
YLIM([LowerYLimit UpperYLimit])
xlabel ('elevation angle (degrees)')
ylabel('actuator displacement (m)')
title('mean')

subplot(2,2,2)
plot(pos,array75, '-o')
YLIM([LowerYLimit UpperYLimit])
xlabel ('elevation angle (degrees)')
ylabel('actuator displacement (m)')
title('75th percentiles')

subplot(2,2,3)
plot(pos, array50, '-o')
YLIM([LowerYLimit UpperYLimit])
xlabel ('elevation angle (degrees)')
ylabel('actuator displacement (m)')
title('50th percentiles')

subplot(2,2,4)
plot(pos, array25, '-o')
YLIM([LowerYLimit UpperYLimit])
xlabel ('elevation angle (degrees)')
ylabel('actuator displacement (m)')
title('25th percentiles')


hold off

% below is where the velocity analysis is conducted.

%here the "velocity" array is initialized, time units will be
%dimensionless, or rather, use the position angle changes as its
%denominator. 

velocityarray = zeros(r,c);
quadratic_derivative_operator=[0 0 0; 2 0 0; 0 1 0];

for j = 1:r
    for i =1:(c-2)
        quadratic_array=zeros(3,3);
        quadratic_array(:,1)=pos(i:1:i+2).^2;
        quadratic_array(:,2)=pos(i:1:i+2);
        quadratic_array(:,3)=[1 1 1]';
        y = act(j,i:1:i+2)';
        quadratic_coefficients = quadratic_array\y;
        velocityvector=quadratic_derivative_operator*quadratic_coefficients;
        if i > 1&&i<(c-2)
            velocityarray(j,i+1)=velocityvector(2)*pos(i+1)+velocityvector(3);
        elseif i==1
            velocityarray(j,i)=velocityvector(2)*pos(i)+velocityvector(3);
            velocityarray(j,i+1)=velocityvector(2)*pos(i+1)+velocityvector(3);
        else
            velocityarray(j,i+1)=velocityvector(2)*pos(i+1)+velocityvector(3);
            velocityarray(j,i+2)=velocityvector(2)*pos(i+2)+velocityvector(3);
        end
    end
end

% The comment below is the function that was used for a linear approximation
% to velocity, if one should wish to bring that back, and a percent difference 
% function to compare linear and quadratic fits. Note that the code boxed below
% has the numbers 2 and 3 after velocityarray, for when both arrays are used.

% velocityarray2 = zeros(r,c);
% for i=1:r
%     for j =1:c
%         if j == 1
%             velocityarray2(i,j)= (act(i,j+1)-act(i,j))/(pos(j+1)-pos(j));
%         elseif j<c
%             velocityarray2(i,j) = (act(i,j+1)-act(i,j-1))/(pos(j+1)-pos(j-1));
%         else
%             velocityarray2(i,j)= (act(i,j)-act(i,j-1))/(pos(j)-pos(j-1));
%         end
%     end
% end
% 
% velocityarray3 = zeros(r,c);
% for i=1:r
%     for j =1:c
%         velocityarray3(i,j)=((velocityarray2(i,j)-velocityarray(i,j))/max([ abs(velocityarray2(i,j)) abs(velocityarray(i,j)) ]))*100;
%     end
% end
% 
% max(max(abs(velocityarray3)))
% mean(mean(abs(velocityarray3)))

figure

maximal_vel_indices=cell(1,c);
minimal_vel_indices=cell(1,c);
maximal_vel_array=zeros(c,c);
minimal_vel_array=zeros(c,c);

%this subplot is for the minima and maxima
for i = 1:c
    [maxi, ndx] = max(velocityarray(:,i));
    maximal_vel_array(:,i)=velocityarray(ndx,:);
    maximal_vel_indices{i}=sprintf('max ar angle %.0f, %.0f',pos(i), ndx);
    [mini, ndy] = min(velocityarray(:,i));
    minimal_vel_array(:,i)=velocityarray(ndy,:);
    minimal_vel_indices{i}=sprintf('min at angle %.0f, %.0f',pos(i), ndy);
end

% the line below is meant to concatenate the matrices so
% that they fit in as a single argument in the legend.
legendarray=[maximal_vel_indices minimal_vel_indices];
plot(pos, maximal_vel_array, '-o', pos, minimal_vel_array, '-o')
legend(legendarray)
xlabel ('elevation angle')
ylabel('position change/angular change')
title(plottitle)

figure

%this subplot is for the average
average = mean(velocityarray);


sortedlist = sort(velocityarray);

%this subplot will be for the 75th percentiles
array75=zeros(c,c);
percentile75=zeros(1,c);
for i = 1:c
    value75=sortedlist(ceil(0.75*r),i);
    k=find(velocityarray(:,i)==value75);
    array75(:,i)=velocityarray(k(1),:);
end

%this subplot will be for the 50th percentiles, it rounds upward
array50=zeros(c,c);
percentile50=zeros(1,c);
for i = 1:c
    value50=sortedlist(ceil(0.50*r),i);
    k=find(velocityarray(:,i)==value50);
    array50(:,i)=velocityarray(k(1),:);
end

%this subplot will be for the 25th percentiles
array25=zeros(c,c);
percentile25=zeros(1,c);
for i = 1:c
    value25=sortedlist(floor(0.25*r),i);
    k=find(velocityarray(:,i)==value25);
    array25(:,i)=velocityarray(k(1),:);
end

UpperYLimit = max([max(max(average)) max(max(array75)) max(max(array50)) max(max(array25))]);
LowerYLimit = min([min(min(average)) min(min(array75)) min(min(array50)) min(min(array25))]);

subplot(2,2,1)
plot(pos, average,'o-')
YLIM([LowerYLimit UpperYLimit])
xlabel ('elevation angle')
ylabel('position change/angular change')
title('mean')

subplot(2,2,2)
plot(pos,array75, '-o')
YLIM([LowerYLimit UpperYLimit])
xlabel ('elevation angle')
ylabel('position change/angular change')
title('75th percentiles')

subplot(2,2,3)
plot(pos, array50, '-o')
YLIM([LowerYLimit UpperYLimit])
xlabel ('elevation angle')
ylabel('position change/angular change')
title('50th percentiles')

subplot(2,2,4)
plot(pos, array25, '-o')
YLIM([LowerYLimit UpperYLimit])
xlabel ('elevation angle')
ylabel('position change/angular change')
title('25th percentiles')

hold off

end
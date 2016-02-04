% File: nrSciGravity.m
%
% Syntax: nrSciGravity(AnalysisCaseDisp,NomActCaseDisp,DatabaseName)
%
% Description:
%       Performs a complete gravity analysis and produces results for a
%       gravitational loading case.  Two input displacements are provided,
%       one for the basis case where the telescope is assumed to be
%       aligned, and a second for the case where the gravity misalignment
%       results are reported.
%
% Input Parameters:
%       AnalysisCaseDisp - the displacement file for the analysis case
%       for which the gravity effects are reported.
%       NomActCaseDisp - the displacement file for the nominal actuator
%       case where the telescope is assumed to be aligned.
%       DatabaseName - the NRCIM database file name
%
% Output Parameters:
%       none.
%
% Required Global Data Structures:
%       none.
%
% Required Data Files:
%       none.
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: nrSciGravity.m,v 1.4 2012/10/30 17:20:32 roberts Exp $
%
% INDENT-OFF*
% $Log: nrSciGravity.m,v $
% Revision 1.4  2012/10/30 17:20:32  roberts
% now prints out command syntax, reordered input arguments
%
% Revision 1.3  2012/10/30 14:55:06  roberts
% removed function argument to exitmatlab on completion
%
% Revision 1.2  2012/10/29 22:38:48  roberts
% now creates files with elevation angle in filename
%
% Revision 1.1  2012/08/02 18:04:04  roberts
% Initial revision
%
% 
%
% INDENT-ON*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%           Herzberg Institute of Astrophysics                  %%%%%
%%%%%%      Astronomy Technology Research Group - Victoria           %%%%%
%
% (c) <2010>				        (c) <2010>
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

function nrSciGravity(AnalysisCaseDisp,NomActCaseDisp,DatabaseName)

% Get the elevation angles from the displacement case file names
ElevBasis = vlUtGetIntFromString(NomActCaseDisp);
ElevPert = vlUtGetIntFromString(AnalysisCaseDisp);

if isempty(ElevBasis) || isempty(ElevPert)
    error('Elevation angle not included in displacement file names');
end

fprintf('Perturbed elevation angle = %d, Basis elevation angle = %d\n',ElevPert,ElevBasis);

% Run nrSciMerit for all cases
fprintf('\n...Calling nrSciMerit(''%s'',''%s'',''%s'')\n',NomActCaseDisp,[num2str(ElevBasis) 'Full'],DatabaseName);
nrSciMerit(NomActCaseDisp,[num2str(ElevBasis) 'Full'],DatabaseName);
fprintf('\n...Calling nrSciMerit(''%s'',''%s'',''%s'')\n',AnalysisCaseDisp,[num2str(ElevPert) 'Full'],DatabaseName);
nrSciMerit(AnalysisCaseDisp,[num2str(ElevPert) 'Full'],DatabaseName);
fprintf('\n...Calling nrSciMerit(''%s'',''%s'',''%s'',''-NomAct'',''%s'')\n',AnalysisCaseDisp,[num2str(ElevPert) 'Align' num2str(ElevBasis)],DatabaseName,[num2str(ElevBasis) 'Full_Data']);
nrSciMerit(AnalysisCaseDisp,[num2str(ElevPert) 'Align' num2str(ElevBasis)],DatabaseName,'-NomAct',[num2str(ElevBasis) 'Full_Data']);

% now generate results
fprintf('\n...Calling vlSciWriteResults(''Gravity'',''%s'',''%s'',''%s'')\n',[num2str(ElevPert) 'Align' num2str(ElevBasis) '_Data.mat'],[num2str(ElevPert) 'Full_Data.mat'],'-PlotAll');
vlSciWriteResults('Gravity',[num2str(ElevPert) 'Align' num2str(ElevBasis) '_Data.mat'],[num2str(ElevPert) 'Full_Data.mat'],'-PlotAll');

% End of nrSciGravity.m

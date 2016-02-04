% File: <nrCimCreateFeaDataScript.m>
%
% Syntax: nrCimCreateFeaDataScript(NrcimInitName, MacroName, FeaDataName)
%
% Description:
%       Writes an Ansys macro that is used to extract CST related nodes and
%       nodal positions from the Ansys model.  Output file name is created
%       from [IM.RunName '_FEADataScript.mac']
%
% Input Parameters:
%       NrcimInitName - Name of .m file to initialize the NRCIM data
%                                   structure.
%       MacroName - Name of macro file to write. Recommend 
%                                   [IM.RunName'_FEADataScript.mac']
%       FeaDataName - Name of FEA Data file to output from Ansys. Recommend
%                                   [IM.RunName'_FEAData']
%
% Output Parameters:
%       None
%
% Required Global Data Structures:
%       NRCIM
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
% $Id: nrCimCreateFeaDataScript.m,v 1.5 2012/07/04 21:03:09 roberts Exp $
%
% INDENT-OFF*
% $Log: nrCimCreateFeaDataScript.m,v $
% Revision 1.5  2012/07/04 21:03:09  roberts
% Updated to include nodal orientation information
%
% Revision 1.4  2010/10/30 18:47:06  roberts
% *** empty log message ***
%
% Revision 1.3  2009/07/08 21:41:31  roberts
% fixed SeqCB surface type
%
% Revision 1.2  2009/07/08 20:48:36  roberts
% updated to new Obj field names in NRCIM
%
% Revision 1.1  2008/04/30 19:36:16  roberts
% Initial revision
%
% INDENT-ON*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%           Herzberg Institute of Astrophysics                  %%%%%
%%%%%%      Astronomy Technology Research Group - Victoria           %%%%%
%
% (c) <2008>				        (c) <2008>
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

function nrCimCreateFeaDataScript(NrcimInit, MacroName, FeaDataName)

global NRCIM

feval(NrcimInit);

fid = fopen(MacroName,'wt');

fprintf(fid,'!! nrCim Ansys Macro to extract nodal interface data from Ansys\n');
fprintf(fid,'!! \tCreated %s\n',date);

%% APDL Initial code

fprintf(fid,'/PREP7\n\n');

fprintf(fid,'!! Set coordinate system to the base cartesian\n');
fprintf(fid,'csys,0\n\n');

fprintf(fid,'!! Write out node numbers and positions\n');
fprintf(fid,'*cfopen,%s,txt	\t! Opens a "command" file\n\n',FeaDataName);

fprintf(fid,'!! Select a node and output to file\n');
fprintf(fid,'allsel,all	    \t! Select everything\n\n');

%% APDL Output nodal data

fprintf(fid,'\n!!!! Starting Node Output    !!!!\n\n');

for SurfIdx = 1:NRCIM.NumSurf
    fprintf(fid,'\n!!!! Surface %d AType = %s ZType = %s\n',SurfIdx,NRCIM.Surf(SurfIdx).AType,NRCIM.Surf(SurfIdx).ZType);
    switch NRCIM.Surf(SurfIdx).ZType
        case 'SeqCB'
            if strncmp(NRCIM.Surf(SurfIdx).ANameSurf,'Name:',5)
                name = NRCIM.Surf(SurfIdx).ANameSurf(6:end);
            else
                error(['Unknown token for sequential surface: ' NRCIM.Surf(SurfIdx).ANameSurf]);
            end
            
            switch NRCIM.Surf(SurfIdx).AType
                case 'Triad'
                    fprintf(fid,'CMSE,S,%s\n',name);
                    fprintf(fid,'CMSE,R,%s\n',NRCIM.Surf(SurfIdx).ANameDisp);
                    fprintf(fid,'CMSE,R,P1\n');
                    local_writenode(fid);
                    fprintf(fid,'CMSE,S,%s\n',name);
                    fprintf(fid,'CMSE,R,%s\n',NRCIM.Surf(SurfIdx).ANameDisp);
                    fprintf(fid,'CMSE,R,P2\n');
                    local_writenode(fid);
                    fprintf(fid,'CMSE,S,%s\n',name);
                    fprintf(fid,'CMSE,R,%s\n',NRCIM.Surf(SurfIdx).ANameDisp);
                    fprintf(fid,'CMSE,R,P3\n');
                    local_writenode(fid);
                case '1Node'
                    fprintf(fid,'CMSE,S,%s\n',name);
                    fprintf(fid,'CMSE,R,%s\n',NRCIM.Surf(SurfIdx).ANameDisp);
                    local_writenode(fid);
                otherwise
                    error('Bad AType');     
            end

        case 'NonSeq'
            if strncmp(NRCIM.Surf(SurfIdx).ANameSurf,'Name:',5)
                nameroot = NRCIM.Surf(SurfIdx).ANameSurf(6:end);
                for ii = 1:NRCIM.Surf(SurfIdx).ZObjNum
                    names{ii} = [nameroot num2str(ii)];
                end
            elseif strncmp(NRCIM.Surf(SurfIdx).AName,'Function:',9)
                funcname = NRCIM.Surf(SurfIdx).AName(10:end);
                names = feval(funcname);
            else
                error('AName token of Name: or Function: was not supplied in NRCIM configuration file');
            end
            
            switch NRCIM.Surf(SurfIdx).AType
                case 'Triad'
                    for ObjIdx = 1:NRCIM.Surf(SurfIdx).ZObjNum
                        fprintf(fid,'CMSE,S,%s\n',names{ObjIdx});
                        fprintf(fid,'CMSE,R,%s\n',NRCIM.Surf(SurfIdx).ANameDisp);
                        fprintf(fid,'CMSE,R,P1\n');
                        local_writenode(fid);
                        fprintf(fid,'CMSE,S,%s\n',names{ObjIdx});
                        fprintf(fid,'CMSE,R,%s\n',NRCIM.Surf(SurfIdx).ANameDisp);
                        fprintf(fid,'CMSE,R,P2\n');
                        local_writenode(fid);
                        fprintf(fid,'CMSE,S,%s\n',names{ObjIdx});
                        fprintf(fid,'CMSE,R,%s\n',NRCIM.Surf(SurfIdx).ANameDisp);
                        fprintf(fid,'CMSE,R,P3\n');
                        local_writenode(fid);
                    end
                case '1Node'
                    for ObjIdx = 1:NRCIM.Surf(SurfIdx).ZObjNum
                        fprintf(fid,'CMSE,S,%s\n',names{ObjIdx});
                        fprintf(fid,'CMSE,R,%s\n',NRCIM.Surf(SurfIdx).ANameDisp);
                        local_writenode(fid);
                    end                    
                otherwise
                    error('Bad AType');  
            end
        otherwise
            error('Bad ZType');
    end
end

%% APDL Ending code

fprintf(fid,'\n\n!!!! Close the file    !!!!\n');
fprintf(fid,'*cfclose\n');

fprintf(fid,'!! Select everything\n');
fprintf(fid,'allsel,all\n');

%% Close the file
fclose(fid);
% End of <nrCimCreateFeaDataScript>

%% function local_writenode

function local_writenode(fid)

%!! Get the node number
fprintf(fid,'nodenum=ndnext(0)\n');

%!! Get the nodal coordinate
fprintf(fid,'nodex = nx(nodenum)\n');
fprintf(fid,'nodey = ny(nodenum)\n');
fprintf(fid,'nodez = nz(nodenum)\n');

%!! Get the nodal orientation
fprintf(fid,'*get,nodethxy,NODE,nodenum,ANG,XY\n');
fprintf(fid,'*get,nodethyz,NODE,nodenum,ANG,YZ\n');
fprintf(fid,'*get,nodethzx,NODE,nodenum,ANG,ZX\n');

fprintf(fid,'*vwrite,nodenum,nodex,nodey,nodez,nodethxy,nodethyz,nodethzx\n');
fprintf(fid,'(F8.0, E16.8, E16.8, E16.8, E16.8, E16.8, E16.8)\n\n');

% end of local_writenode
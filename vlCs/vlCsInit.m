%   File: vlCsInit.m 
%
%   Syntax: vlCsInit
%
%   Discussion:
%       Fills CST global data structure.  All changes to the CST data
%       structure should be made here.
%       
%   Input Parameters:
%	zgTransform - Filename of ANSYS CSYS0 to Z transform
%      distFilename - Filename of model specific distances function
%
%   Output Parameters:
%       Modifies global data structure CST
%
% Required Global Data Structures:
%       CST
%	OC 
%	SD 
%	IM 
%
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
% static char rcsid[] = "$Id: vlCsInit.m,v 1.12 2008/02/12 22:53:56 msmith Exp $";
% INDENT-OFF*
% $Log: vlCsInit.m,v $
% Revision 1.12  2008/02/12 22:53:56  msmith
% Added code to handle distFilename being passed in as a parameter.
% Defaults to vlCsM1pGetDistances_3mirr for backwards compatibility.
%
% Revision 1.11  2006/11/23 03:12:27  msmith
% Added distFilename parameter to allow for specification of distances file.
%
% Revision 1.10  2006/11/08 22:19:23  msmith
% Changed "filename" parameter to "zgTransform".  Updated parameter comment.
%
% Revision 1.9  2006/01/05 20:34:06  msmith
% Restructured code to make SD.FEAInterface query simpler.
%
% Revision 1.8  2004/07/28 20:48:56  kerleyd
% modified call to vlCsHexMirr
%
% Revision 1.7  2004/07/13 18:42:01  kerleyd
% modified to use new SD structure to create the coordinate transforms
% in the in 'ring' order
%
% Revision 1.6  2004/06/08 22:03:37  msmith
% Changed OC.PrimRad to more descriptive OC.PrimRadCurv.
%
% Revision 1.5  2003/11/24 22:24:31  stretchn
% changed 'string' to 'num2str'
%
% Revision 1.4  2003/10/28 22:43:08  stretchn
% Now you pass in the filename, instead of accessing IM
%
% Revision 1.3  2003/10/27 20:08:24  stretchn
% Added 'tic' command to make time printout correct.
%
% Revision 1.2  2003/09/18 23:43:49  stretchn
% Updated header
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

function vlCsInit(zgTransform, distFilename)

global CST
global OC
global SD

disp('Entering vlCsInit ...');

tic;

%
% Store distances in CST.
% Legacy code did not specify distFilename or use CST.Distances.
% Use vlCsM1pGetDistances_3mirr as the default (MF Rel 3-5).
%

if nargin < 2
    distFilename = 'vlCsM1pGetDistances_3mirr';
end

CST.Distances = eval(distFilename);


%
% Z to G coordinate transformation
%

CST.ZG = eval(zgTransform);


%
% Defines x, y and z columns in the interface file
%

xyzColumn = [vlSdGetInterColumn('xPos'),vlSdGetInterColumn('yPos'),...
                vlSdGetInterColumn('zPos')];
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sequential CS     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Sequential Mirror Surface Coordinate System
%

for i = 1:OC.NumSeqs
% strrep(source, search, rep)
% OC.SeqSurfName{i} - 'Secondary', 'Tertiary', 'Focus'
% OC.ZMSeqCSTName{i} - 'vlCsRcZsecondary_3mirr', 'vlCsRcZtertiary_3mirr', 'vlCsRcZfocus_3mirr'
% convert name to lower case then substitute 'GetDistances' with 'Z' + name

    CST.ZMSeq(:,:,i) = feval(OC.ZMSeqCSTName{i});
end

%
% Sequential Mirror Cell Coordinate System
%

optSeq = SD.SSoutput(SD.SeqMap(1:3:end));

nodeColumnIndex = vlSdGetInterColumn('nodeNum');
indexColumnIndex = vlSdGetInterColumn('index');

for surf = 1:OC.NumSeqs
    
    %
    % Gets traid of each sequential surface 
    %
    
    for triad = 1:3
        node(triad) = SD.FEAInterface(SD.FEAInterface(:,nodeColumnIndex) ...
                   == optSeq((surf-1)*3+triad), indexColumnIndex);
    end 
    
    %    
    % Gets coordinates for each triad 
    %
    
    P1 = vlCsPMult(CST.ZG,SD.FEAInterface(node(1),xyzColumn)');
    P2 = vlCsPMult(CST.ZG,SD.FEAInterface(node(2),xyzColumn)');
    P3 = vlCsPMult(CST.ZG,SD.FEAInterface(node(3),xyzColumn)');
    
    CST.ZCSeq(:,:,surf) = vlCsTriadTrans(P1,P2,P3);
end

%
% Sequential Cell to Mirror CST (CST.CMSeq)
%

for i = 1:OC.NumSeqs
    CST.CMSeq(:,:,i) = vlCsMult(vlCsInv(CST.ZCSeq(:,:,i)),CST.ZMSeq(:,:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-Sequential CS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Primary Mirror Surface coordinate System (CST.ZMNonSeq)
%

primSegId = vlSdGetIdentNum('PrimarySeg');
segSupport = vlSdGetSupportNum(primSegId,'segNum');
ringSupport = vlSdGetSupportNum(primSegId,'ringNum');
existsSupport = vlSdGetSupportNum(primSegId,'exists');

for segNum = 1:OC.NumNonSeqs
    
    %
    % checks if segment is missing, if missing segment fill in dummy data
    %
    
    if SD.FEAInterface(SD.FEAInterface(:,...
            vlSdGetInterColumn('id'))==primSegId & SD.FEAInterface(:,...
            vlSdGetInterColumn(['S' num2str(segSupport)]))==segNum,...
            vlSdGetInterColumn(['S' num2str(existsSupport)]))
        
        %
        % returns ring # for a give seg #
        %
        
        ringNum = SD.FEAInterface(SD.FEAInterface(:,...
                vlSdGetInterColumn('id'))==primSegId & SD.FEAInterface(:,...
                vlSdGetInterColumn(['S' num2str(segSupport)]))==segNum,...
                vlSdGetInterColumn(['S' num2str(ringSupport)]));
        
        [T,ETseg] = vlCsHexMirr(segNum,ringNum,OC.SegSize,OC.PrimCC,OC.PrimRadCurv);
        CST.ZMNonSeq(:,:,segNum) = T;
        
    else
        
        %
        % If missing segment,identity matrix transform is used 
        %
        
        CST.ZMNonSeq(:,:,segNum) = eye(3,4);
    end
end    

%
% Primary Mirror cell coordinate system  (CST.ZCNonSeq)
%

optNonSeq = SD.SSoutput(SD.NonSeqMap(1:3:end));

for segNum = 1:OC.NumNonSeqs
	
    %
    % checks if segment is missing, if missing segment fill in dummy data
    %
    
    if SD.FEAInterface(SD.FEAInterface(:,...
        vlSdGetInterColumn('id'))==primSegId & SD.FEAInterface(:,...
        vlSdGetInterColumn(['S' num2str(segSupport)]))==segNum,...
        vlSdGetInterColumn(['S' num2str(existsSupport)]))

        %   
        % Gets traid of each non-sequential surface 
        %

        for triad = 1:3

            %
            % return interface file index for a given node #
            %

            node(triad) = SD.FEAInterface(SD.FEAInterface(:,...
                vlSdGetInterColumn('nodeNum'))==optNonSeq((segNum-1)*3+triad),...
                vlSdGetInterColumn('index'));
        end

        %
        % Gets coordinates for each triad 
        % use the triad points directly to calculate the transform...
        %

        P1 = vlCsPMult(CST.ZG,SD.FEAInterface(node(1),xyzColumn)');
        P2 = vlCsPMult(CST.ZG,SD.FEAInterface(node(2),xyzColumn)');
        P3 = vlCsPMult(CST.ZG,SD.FEAInterface(node(3),xyzColumn)');        

        CST.ZCNonSeq(:,:,segNum) = vlCsTriadTrans(P1,P2,P3); 

      else  
          
        %
        % If missing segment,identity matrix transform is used 
        %
          
        CST.ZCNonSeq(:,:,segNum) = eye(3,4);         
        
      end
end

%
% Primary Cell to Mirror CST (CST.CMNonSeq)
%

for i = 1:OC.NumNonSeqs
    CST.CMNonSeq(:,:,i) = vlCsMult(vlCsInv(CST.ZCNonSeq(:,:,i)),CST.ZMNonSeq(:,:,i));
end

disp(['    ... Done, elapsed time = ' num2str(toc)]);

% End of vlCsInit

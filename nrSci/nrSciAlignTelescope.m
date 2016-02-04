% File: <nrSciAlignTelescope.m>
%
% Syntax: nrSciAlignTelescope
%
% Description:
%       Routine to find fit of displaced telescope nodes to align the
%       primary mirror
%
% Input Parameters:
%       none
%
% Output Parameters:
%       none.  (Modifies the RES data structure)
%
% Required Global Data Structures:
%         global CST
%         global RES
%         global NRCIM
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
% $Id: nrSciAlignTelescope.m,v 1.8 2012/07/24 19:56:37 roberts Exp $
%
% INDENT-OFF*
% $Log: nrSciAlignTelescope.m,v $
% Revision 1.8  2012/07/24 19:56:37  roberts
% updated DebugMode behavior
%
% Revision 1.7  2011/03/08 22:32:01  roberts
% fixed header
%
%
% Revision 1.6. Added code to compare average segment position with LS solution.  Removed
% 100x angle scaling for fminsearch.  Changed tolerances on fminsearch.
%
% INDENT-ON*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%           Herzberg Institute of Astrophysics                  %%%%%
%%%%%%      Astronomy Technology Research Group - Victoria           %%%%%
%
% (c) <2009>				        (c) <2009>
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

function nrSciAlignTelescope

% Create "module" variables that are shared among the nested functions.
A = [];
B = [];
X = [];

doAlignment;


%% doAlignment : Main function

function doAlignment

    % declare global variables
    global CST
    global RES
    global NRCIM
    global OC
    global MA

    %% Determine the method used to align the telescope
    % 1 = LS   = Least squares calculation of transform T
    % 2 = LS+  = Minimum actuator stroke
    % 3 = Keck = Keck planar method
    % 4 = AvXY = Uses the average XY of the
    %                segments then calculates minimum
    %                actuator stroke using remaining DOFs

    if isfield(RES, 'ActFit') && isfield(RES, 'ActFitNdx')
      fitMethod = RES.ActFitNdx;  % Assume that it is set correctly.
    else
      warning('RES.ActFit or RES.ActFitNdx not set.  Using LS+ method.');
      fitMethod = 2;  % Set the default method here.
    end
    
    %% Create two vectors of actuator positions in the Z coordinate system
    % A is a vector of nominal actuator positions
    % B is a vector of displaced actuator positions

    M1_surfnum = 1;
    NumSegs = NRCIM.SegsPerSurf(M1_surfnum);
    
    % no longer have SD data structure, so calculate value used for
    % SD.TriadSizeNonSeq
    
    M1TriadSizes = [CST.Meta([CST.Meta(:).NrcimSurf]==1).ADim];
    M1TriadSize = mean(M1TriadSizes);
    M1TriadSizeTol = (max(M1TriadSizes) - min(M1TriadSizes))/M1TriadSize;
    if M1TriadSizeTol > 1e-6
        warning('Triad sizes for M1 segments vary more than expected, tol = %f\n',M1TriadSizeTol);
    end
    
    % The nominal positions of the actuators in the C coordinate system are
    %             [-sqrt(3)/2*M1TriadSize -M1TriadSize/2 0]';
    %             [ sqrt(3)/2*M1TriadSize -M1TriadSize/2 0]';
    %             [ 0                             M1TriadSize   0]';
    
    % Define the nominal actuator positions in the C coordinate system
    TSizeHalf = M1TriadSize / 2;
    TSizeHalfR3 = TSizeHalf * sqrt(3);
    
    % Read variables as C coordinate, nominal P1 (P1 == first actuator)
    C_NP1 = [-TSizeHalfR3 -TSizeHalf 0]';
    C_NP2 = [ TSizeHalfR3 -TSizeHalf 0]';
    C_NP3 = [ 0 M1TriadSize 0]';
    
    A = zeros(NumSegs*9,1);
    B = zeros(NumSegs*9,1);
    
    % Define an index to fill A and B
    Aidx = 0;
    
    for ii = 1:CST.Num
        % Does the CST belong to the primary mirror?
        if CST.Meta(ii).NrcimSurf == M1_surfnum
            
            % Read the displaced actuator positions in nodal coordinates.
            % Note that in Ansys the nodal coordinates are aligned with the
            % C CST for each segment. 
            C_DP1 = RES.Y(CST.Meta(ii).YMap(1:3));
            C_DP2 = RES.Y(CST.Meta(ii).YMap(4:6));
            C_DP3 = RES.Y(CST.Meta(ii).YMap(7:9));
            
            % Compute the nominal actuator positions in the Z CST for each
            % segment and store then in A
            A(Aidx+1:Aidx+3) = vlCsPMult(CST.ZC(:,:,ii),C_NP1);
            A(Aidx+4:Aidx+6) = vlCsPMult(CST.ZC(:,:,ii),C_NP2);
            A(Aidx+7:Aidx+9) = vlCsPMult(CST.ZC(:,:,ii),C_NP3); 
             
            % Compute the displaced actuator positions in the Z CST for
            % each segment and store them in B
            B(Aidx+1:Aidx+3) = vlCsPMult(CST.ZC(:,:,ii),C_NP1+C_DP1);
            B(Aidx+4:Aidx+6) = vlCsPMult(CST.ZC(:,:,ii),C_NP2+C_DP2);
            B(Aidx+7:Aidx+9) = vlCsPMult(CST.ZC(:,:,ii),C_NP3+C_DP3);  
                                        
            % Increment index
            Aidx = Aidx + 9;
        end
    end
        

    % If a nominal actuator position vector was supplied using the
    % '-NomAct' flag, then use it instead.
    if ~isempty(RES.NomActPosFile)
        A = RES.ActPosNom.AlignZ;
    end
    
    % save the A vector used above
    RES.ActPosNom.A = A;

    % Save the displaced actuator positions
    RES.ActPosDisp = B;
    
    if RES.DebugMode
        Alist = reshape(A,3,length(A)/3)';
        Blist = reshape(B,3,length(B)/3)';
        Dlist = Blist - Alist;
        save 'PrimOrig.txt' Alist -ASCII;
        save 'PrimDisp.txt' Dlist -ASCII;
    end    


    if fitMethod == 3 % Keck

       [ coefs, euler, T, rmsErr ] = local_Keck;  % what is flag?

       fprintf('Keck z=a+bx+cy:  a=%.4e b=%.4e c=%.4e\n',coefs(1),coefs(2),coefs(3));
       fprintf('Keck Fitting error: %.4e\n', rmsErr);

       RES.EULZZP = euler;
       RES.CSTZZP = T;  % vlCsInvEulerXYZ(euler);
       RES.CSTZPZ = vlCsInv(RES.CSTZZP);

       % Done least squares fitting
       fprintf('....done in %.1f Seconds\n\n',toc);

    %%  Least Squares Fit to Primary, Calc. of Actuator and center of curvature
    
    elseif ismember(fitMethod, [ 1, 2]) % All others

        fprintf('Calculating alignment of primary mirror\n');

        % Now we want to find a CST that best fits the data 
        
        % save everything for debugging
        if RES.DebugMode
            save(strcat(RES.outFileRoot,'_AllData'));
        end
        
        % The x vector consists of X,Y,Z translations followed by Rx,Ry,Rz
        % rotatons. Angles are in degrees
        % In order to have the same magnitude on the input values, we can
        % conider scaling angles with respect to the displacements.
        % In degrees, the angles are 4 times less sensitive than
        % displacements.  I.E. sin(4deg)*15 m = 1m. This is close enough to
        % not scale.
        

        if ismember(fitMethod, [ 1, 2]) % LS or LS+
            % As a starting point find the least squares solution
            [x,fval,exitflag,output] = fminsearch(@local_lscalc,[0 0 0 0 0 0], ...
            optimset('TolX',1e-8,'TolFun',1e-7,'MaxFunEvals',4000,'MaxIter',4000));

            % Set the module variable X equal to the solution from the 
            % least squares fit
            if exitflag ~= 1 
                error('fminsearch did not converge');
            end
            fprintf('Least squares solution complete\nx = %.4e %.4e %.4e %.4e %.4e %.4e fval = %.4e\n',x(1),x(2),x(3),x(4),x(5),x(6),fval);
        end
        
        %% This is code to compare the average segment position with the LS
        % solution
        
        if ismember(fitMethod, [1, 2]) % LS or LS+ % AvXY
            % Use the average position of the segments for the XY position
            
            % For all the actuator positions in B calculate the segment
            % position and average the XYZ
            
            T_ZZP = vlCsInvEulerXYZ(x);
            
            Idx = 0;
            Ptotal = zeros(3,1);
            
            for ii=1:NumSegs
                P1 = B(Idx+1:Idx+3);
                P2 = B(Idx+4:Idx+6);
                P3 = B(Idx+7:Idx+9);
                
                % Calculate the transform to the triad center CST ZC
                T_ZC = vlCsTriadTrans(P1,P2,P3);
                T_ZPC = vlCsMult(vlCsInv(T_ZZP),T_ZC);
                
                % Calculate the position of the segment vertex CST ZC x Point in C
                % Recall that Z axis of Z and C point away from sky, so
                % PrimOffset is negative
                
                P = vlCsPMult(T_ZPC,[0 0 -OC.PrimOffset]');
                Ptotal = Ptotal+P;
                
                Idx = Idx+9;
            end
            % Calculate average position
            Ptotal = Ptotal / NumSegs;
            
            % Set X positions
            fprintf('Differential of Ave XY to LS solution = %f,%f\n',Ptotal(1),Ptotal(2));
        end

        %% Perform actuator stroke minimization if required.
        if fitMethod == 2 % LS+

            % Now minimize the actuator displacement rather than 
            % least squares for minimum.
            % First optimize for only Rx and Z, then for Rx,Ry and Z 
            % actuator stroke but only with Rx, Ry and Z as variables

            X = x;
            [x,fval,exitflag,output] = ...
            fminsearch(@local_minzcalc_zrx,[x(3) x(4)],...
            optimset('TolX',1e-8,'TolFun',1e-7,'MaxFunEvals',4000,'MaxIter',4000));

            if exitflag ~= 1 
                error('fminsearch did not converge');
            end

            % reconstitute x
            x = [X(1) X(2) x(1) x(2) X(5) X(6)];

            fprintf('Min Rx,Z solution complete\nx = %.4e %.4e %.4e %.4e %.4e %.4e  fval = %.4e\n',x(1),x(2),x(3),x(4),x(5),x(6),fval);

            X = x;
            [x,fval,exitflag,output] = ...
            fminsearch(@local_minzcalc_zrxy,[x(3) x(4) x(5)],...
            optimset('TolX',1e-8,'TolFun',1e-7,'MaxFunEvals',4000,'MaxIter',4000));
            if exitflag ~= 1 
                error('fminsearch did not converge');
            end

            % reconstitute x
            x = [X(1) X(2) x(1) x(2) x(3) X(6)];

            fprintf('Min Rx,Ry, Z solution complete\nx = %.4e %.4e %.4e %.4e %.4e %.4e  fval = %.4e\n',x(1),x(2),x(3),x(4),x(5),x(6),fval);

        end % fitmethod 2

        RES.EULZZP = x;
        RES.CSTZZP = vlCsInvEulerXYZ(x);
        RES.CSTZPZ = vlCsInv(RES.CSTZZP);

        % Done actuator fitting
        fprintf('....done in %.1f Seconds\n\n',toc);

    end  % fitMethod=1,2  (check that this is correct location)

%% Process the output data
% The following outputs are calculated:
% 
% The outputs needed are:
% •	A list displaced actuator positions in the ZP coordinate system
% [ActPosAlignZ]
% •	The displaced actuator positions in the C coordinate systems relative
% to ZP [ActPosAlignC].  
% •	The displaced actuator positions in the Nodal coordinate systems This
% is equivalent to the residual actuator stroke that is required. (or is
% it?)  [ActPosAlignN] 
% •	A Y vector in nodal coordinates based on the new ZP coordinate system
% [ActPosAlignY]
% •	The actuator stroke required to align the primary mirror.  This is the
% Z-value from each actuator nodal displacement. [ActStroke]
% 
% Note: SSOutput is in the order of the Y vector

% Initialize data structures
ActPosAlignZ = zeros(NumSegs*9,1);
ActPosAlignC = zeros(NumSegs*9,1);
ActPosAlignN = zeros(NumSegs*9,1);
ActPosAlignY = zeros(CST.Num*9,1);

ActIdx = 0;

for ii = 1:CST.Num
    % Does the CST belong to the primary mirror?
    if CST.Meta(ii).NrcimSurf == M1_surfnum
        ActPosAlignZ(ActIdx+1:ActIdx+3) = vlCsPMult(RES.CSTZPZ,RES.ActPosDisp(ActIdx+1:ActIdx+3));
        ActPosAlignZ(ActIdx+4:ActIdx+6) = vlCsPMult(RES.CSTZPZ,RES.ActPosDisp(ActIdx+4:ActIdx+6));
        ActPosAlignZ(ActIdx+7:ActIdx+9) = vlCsPMult(RES.CSTZPZ,RES.ActPosDisp(ActIdx+7:ActIdx+9));
        
        ActPosAlignC(ActIdx+1:ActIdx+3) = vlCsPMult(vlCsInv(CST.ZC(:,:,ii)),ActPosAlignZ(ActIdx+1:ActIdx+3));
        ActPosAlignC(ActIdx+4:ActIdx+6) = vlCsPMult(vlCsInv(CST.ZC(:,:,ii)),ActPosAlignZ(ActIdx+4:ActIdx+6));
        ActPosAlignC(ActIdx+7:ActIdx+9) = vlCsPMult(vlCsInv(CST.ZC(:,:,ii)),ActPosAlignZ(ActIdx+7:ActIdx+9));
        
        ActPosAlignN(ActIdx+1:ActIdx+9) = ActPosAlignC(ActIdx+1:ActIdx+9) - [C_NP1 ; C_NP2 ; C_NP3];

        ActPosAlignY(CST.Meta(ii).YMap(1:9)) = ActPosAlignN(ActIdx+1:ActIdx+9);
        
        ActIdx = ActIdx+9;
    end
end

    RES.ActPosNom.AlignZ = ActPosAlignZ;
        
    ActStroke = zeros(NumSegs*3,1);
    for ii = 1:NumSegs*3
        Idx = (ii-1)*3+1;
        ActStroke(ii) = vlSciActStroke(RES.ActPosNom.A(Idx:Idx+2),ActPosAlignZ(Idx:Idx+2));
    end
    
    RES.ActPosAlignC = ActPosAlignC;
    RES.ActPosAlignN = ActPosAlignN;
    RES.ActPosAlignY = ActPosAlignY;
    RES.ActStroke = ActStroke;

    % Done transforming points
    fprintf('....done in %.1f Seconds\n\n',toc);
end % doAlignment


%%      local_lscalc

function res = local_lscalc(eulerSearchPosition)
    
    T = vlCsInvEulerXYZ(eulerSearchPosition);

    % transform A with T
    AP = zeros(size(A));

    for ii = 1:3:length(A)
        AP(ii:ii+2) = vlCsPMult(T,A(ii:ii+2));
    end

    % least squares solution...
    res = sum((AP-B).^2);

end   % local_lscalc

%%      local_minzcalc_zrxy

function res = local_minzcalc_zrxy(inp)

    % X,Y and RZ are stored in the module X structure
    %inp(1) is Z
    %inp(2) is RX
    %inp(3) is RY
    
    angleset = [X(1) X(2) inp(1) inp(2) inp(3) X(6)];
    
    T = vlCsInvEulerXYZ(angleset);

    % transform A with T

    AP = zeros(size(A));
    actstroke = zeros(size(A));

    for ii = 1:3:length(A)
        AP(ii:ii+2) = vlCsPMult(T,A(ii:ii+2));
        actstroke(ii) = vlSciActStroke(AP(ii:ii+2),B(ii:ii+2));
    end
    res = max(abs(actstroke));
     
    %fprintf('res = %.4e\n',res);
    
end   % local_minzcalc_zrxy

%%      local_minzcalc_zrxy

function res = local_minzcalc_zrx(inp)

    % X,Y and RZ are stored in the module X structure
    %inp(1) is Z
    %inp(2) is RX
    
    angleset = [X(1) X(2) inp(1) inp(2) X(5) X(6)];
    
    T = vlCsInvEulerXYZ(angleset);

    % transform A with T

    AP = zeros(size(A));
    actstroke = zeros(size(A));

    for ii = 1:3:length(A)
        AP(ii:ii+2) = vlCsPMult(T,A(ii:ii+2));
        actstroke(ii) = vlSciActStroke(AP(ii:ii+2),B(ii:ii+2));
    end
    res = max(abs(actstroke));
     
    %fprintf('res = %.4e\n',res);
    
end   % local_minzcalc_zrx


%% local_Keck

function [ coefs, euler, T, rmsErr ] = local_Keck  % what is flag?

    % A and B are column vectors so we extract the x and y elements as
    % columns of M
    
    fprintf('Computing Keck attitude plane\n');
    
    xCol = B(1:3:end);
    yCol = B(2:3:end);
    zCol = B(3:3:end) - A(3:3:end);  % Should this be B-A or A-B?
    n = numel(xCol);

    M = [ ones(n,1) xCol yCol ];
    Minv = pinv(M);
    
    % Compute [ a b c ] of the Keck note
    coefs = Minv * zCol;
    a = coefs(1);
    b = coefs(2);
    c = coefs(3);

    
    % Compute RMS error
    fitted = M * coefs;
    rmsErr = norm(fitted-zCol);
    
    if RES.DebugMode
        fprintf('rmsErr(fitted-zCol) = %f\n', rmsErr);
        fprintf('rmsErr(A+fitted-B) = %f\n', norm(A(3:3:end)+fitted-B(3:3:end)));
    end

    % Convert slopes and offset into a transform.
    % A positive rotation about the X axis corresponds to a positive Z.
    % A negative rotation about the Y axis corresponds to a positive Z.
    rx = atan(c) * 180 / pi;
    ry = -atan(b) * 180 / pi;
    euler = [ 0 0 a rx ry 0 ];  % need to check if the signs are correct for these
    T = vlCsInvEulerXYZ(euler);
    
end   % local_Keck

end  % nrSciAlignTelescope

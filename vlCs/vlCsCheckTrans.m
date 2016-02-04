% File: vlCsCheckTrans.m
%
% Syntax: [ status, xyzSeg, xyzTriad ] = vlCsCheckTrans(model, ...
%                                                       chkZemax, doGraphs)
%
% Description:
%       The vlCsCheckTrans function performs a variety of checks on the
%       transforms stored in the CST structure.  These checks include the
%       following.
%
%    1. Compare actual and expected dimensions of the transforms.
%    2. Determine if ZEMAX prescription uses segments with equal sized 
%       projections.
%    3. Check the ZEMAX transforms against the IM ZMSeq transforms.
%    4. Compare the CST.ZMSeq transforms with the CST.ZCSeq transforms.
%    5. Compare the CST.ZMNonSeq with the CST.ZCNonSeq.
%    6. Compare the CST.ZMNonSeq transforms with ZEMAX segment transforms.
%    7. Compute the maximum distance between segment centers (from SD struct),
%       mirror centers (from CST.ZMNonSeq transforms) and triad centers 
%       (from SD structure).
%
% Input Parameters:
%       model    - (string) Short "name" for current model (e.g. "M1+")
%       chkZemax - (scalar) Validate against ZEMAX if chkZemax is nonzero.
%       doGraphs - (scalar) Controls display and saving of graphs.
%                   0 - Do not display or save graphs.
%                   1 - Display but do not save graphs.
%                   2 - Display and save graphs.
%
% Output Parameters:
%       status   - (Data Dimension) Description
%                   continues here.
%       xyzSeg   - (OC.NumNonSeqs * 4) Position (in M1CS) of each segment.
%                  Column 1..3 = X, Y and Z positions respectively
%                  Column 4 = 1 if segment exists, 0 if segment is missing.
%       xyzTriad - (3*N, 5) Position and identification of triad points.
%                  (Three rows for each non-missing segment.)
%                  Column 1..3 = X, Y and Z positions respectively
%                  Column 4 = Segment number of triad point.
%                  Column 5 = Number of triad point within triad (segment).
%
% Required Global Data Structures:
%       CST - This contains the transforms we are checking.
%       OC  - Contains number of mirrors, missing segments, etc.
%       SD  - Contains triad and segment positions, etc.
%
% Required Data Files:
%       N/A
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: vlCsCheckTrans.m,v 1.17 2008/02/15 18:42:01 msmith Exp $
%
% INDENT-OFF*
% $Log: vlCsCheckTrans.m,v $
% Revision 1.17  2008/02/15 18:42:01  msmith
% Updated checking of segment C-M distances for MF 7.0.
%
% Revision 1.16  2008/02/09 01:12:52  msmith
% Improved error message on C-M distances too small/large detection.
%
% Revision 1.15  2006/12/07 22:03:16  msmith
% Improved error checking of CST.Distances = IM.GetDistancesFile validation.
%
% Revision 1.14  2006/11/23 03:13:28  msmith
% Added code to check if CST.Distances agrees with the values specified in
% IM.GetDistancesFile (if they both exist).
%
% Revision 1.13  2006/05/23 20:04:27  msmith
% Removed commented code which compares C-M distance to SD.TriadSizeSeq.
% Code would not work under new SD.TriadSizeSeq and scaling factors were
% too arbitrary to be useful.
%
% Revision 1.12  2006/03/20 21:23:15  msmith
% Corrected copy/paste error in C-M distances too large message.
%
% Revision 1.11  2006/03/17 19:04:42  msmith
% Updated validation tests: refined CS alignment tests, adjusted tolerances.
%
% Revision 1.10  2005/03/24 00:07:17  msmith
% Moved the assignment to con and seg above the if Zemax test so that the
% variables are always defined.
%
% Revision 1.9  2005/02/16 23:34:15  msmith
% Changed zmxSeqThreshold from 1e-9 to 1e-10.
% Added printing of error value in ZMSeq = ZCSeq * CMSeq check.
% Added code to check M-C distances for all existing primary mirror segments.
% Added code to store maximum distances between segment centers in status
%   result rather than just displaying it.
%
% Revision 1.8  2005/02/15 00:10:34  msmith
% Added zmxSeqThreshold and zmxNonSeqThreshold variables.
% Added code to check sizes of transformation arrays.
% Improved status messages (include error values, limit # of errors).
% Added code to detect anti-parallel C and M for sequentials.
% Changed zemaxNorm to zemaxErr and scaled zemaxErr.
% Added missing "abs" calls on tests of Z for opposing segments.
%
% Revision 1.7  2005/02/04 20:16:08  msmith
% Fixed index in Ck and Mk messages (now C{k+1} and M{k+1}).
% Changed threshold for CST.CMSeq(3,3,k) from 0.99 to 0.999.
% Added test of Ck to Mk distance based on CST.CMSeq(3,4,k).
% Changed threshold for CST.CMNonSeq calculation from 1.0e-2 to 1.0e-15.
% Added printing of maximum error from CST.CMNonSeq check.
%
% Revision 1.6  2005/02/03 22:04:22  msmith
% Added sanity check on CST.CMSeq.
% Added warning messages to returned list of status messages.
% Added sanity check on CST.CMNonSeq.
% Added code to test alignment of C and M using CMSeq.
%
% Revision 1.5  2005/01/28 00:52:17  msmith
% Replaced st_* constants with pfx and call to vlUtValPrefixes.
% Added/updated comments in header.
%
% Revision 1.4  2004/12/01 19:05:48  msmith
% Added checks to compare CST.ZMNonSeq with ZEMAX segment transforms.
%
% Revision 1.3  2004/11/19 00:31:39  msmith
% Fixed equal projection area test (added missing diff call).
% Changed cell array from two columns to one column.
%
% Revision 1.2  2004/11/11 00:22:36  msmith
% Initial completed transform checking routine.
%
% Revision 1.1  2004/10/04 20:33:57  dunn
% Initial revision
%
% INDENT-ON*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%           Herzberg Institute of Astrophysics                  %%%%%
%%%%%%      Astronomy Technology Research Group - Victoria           %%%%%
%
% (c) <2003>                                  (c) <2003>
% National Research Council              Conseil national de recherches
% Ottawa, Canada, K1A 0R6                  Ottawa, Canada, K1A 0R6
% All rights reserved                        Tous droits reserves
%                                         
% NRC disclaims any warranties,         Le CNRC denie toute garantie
% expressed, implied, or statu-         enoncee, implicite ou legale,
% tory, of any kind with respect        de quelque nature que se soit,
% to the software, including            concernant le logiciel, y com-
% without limitation any war-           pris sans restriction toute
% ranty of merchantability or           garantie de valeur marchande
% fitness for a particular pur-         ou de pertinence pour un usage
% pose.  NRC shall not be liable        particulier.  Le CNRC ne
% in any event for any damages,         pourra en aucun cas etre tenu
% whether direct or indirect,           responsable de tout dommage,
% special or general, consequen-        direct ou indirect, particul-
% tial or incidental, arising           ier ou general, accessoire ou
% from the use of the software.         fortuit, resultant de l'utili-
%                                       sation du logiciel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ status, xyzSeg, xyzTriad ] = vlCsCheckTrans(model, chkZemax, doGraphs)

global CST
global OC
global SD
global IM

status = [];

pfx = vlUtValPrefixes;

testStatus = cell(1,1);

zmxSeqThreshold = 1.0e-10;
zmxNonSeqThreshold = 1.0e-9;


if chkZemax == 0

    msg = 'Not performing ZEMAX tests (chkZemax=0).';
    testStatus{1,1} = [ pfx.info msg ];
    status = [ status ; testStatus ];

else

    if ~ispc
        msg = 'ZEMAX only supported on Windows.  Zemax tests are disabled.';

        warning(msg);
        chkZemax = 0;

        testStatus{1,1} = [ pfx.info msg ];
        status = [ status ; testStatus ];

    else

        % Open a ZEMAX DDE connection if required.

        zemaxStatus = vlMzConnected;

        if (zemaxStatus ~= 1)
            vlMzOpen;
        end

        % Ensure that ZEMAX DDE connection exists.

        zemaxStatus = vlMzConnected;

        if (zemaxStatus ~= 1)

            msg = 'Could not connect to ZEMAX.  Zemax tests are disabled.';

            warning(msg);
            chkZemax = 0;

            testStatus{1,1} = [ pfx.info msg ];
            status = [ status ; testStatus ];

        end
    end
end


% Check the sizes of the CST arrays.

st = checkTransformSize('CST.ZG', size(CST.ZG), [ 3, 4 ]);
status = [ status ; st ];

st = checkTransformSize('CST.ZMSeq', size(CST.ZMSeq), [ 3, 4, OC.NumSeqs ]);
status = [ status ; st ];

st = checkTransformSize('CST.ZCSeq', size(CST.ZCSeq), [ 3, 4, OC.NumSeqs ]);
status = [ status ; st ];

st = checkTransformSize('CST.CMSeq', size(CST.CMSeq), [ 3, 4, OC.NumSeqs ]);
status = [ status ; st ];

st = checkTransformSize('CST.ZMNonSeq', size(CST.ZMNonSeq), ...
                        [ 3, 4, OC.NumNonSeqs ]);
status = [ status ; st ];

st = checkTransformSize('CST.ZCNonSeq', size(CST.ZCNonSeq), ...
                        [ 3, 4, OC.NumNonSeqs ]);
status = [ status ; st ];

st = checkTransformSize('CST.CMNonSeq', size(CST.CMNonSeq), ...
                        [ 3, 4, OC.NumNonSeqs ]);
status = [ status ; st ];

% Check the CST.Distances fields against IM.GetDistancesFile
% A warning is not issued for models prior to version 6.0 since
% the GetDistancesFile field did not exist.


if vlSciMeritGetVersion >= 6.0

    if ~isfield(IM, 'GetDistancesFile')

	warning('IM structure missing field GetDistancesFile');
	testStatus{1,1} = [ pfx.fail msg ];
	status = [ status ; st ];

    elseif ~isfield(CST, 'Distances')

	warning('CST structure missing field Distances');
	testStatus{1,1} = [ pfx.fail msg ];
	status = [ status ; st ];

    else
        % Now check that the values in CST.Distances match those returned
        % by the function specified in IM.GetDistancesFile.

	D = eval(IM.GetDistancesFile);

	Fnames = fieldnames(D);

	for i=1:length(Fnames)
	    d1 = CST.Distances.(Fnames{i});
	    d2 = D.(Fnames{i});
	    if d1 ~= d2
		msg = 'CST.Distances.%s = %.6g != %.6g';
		warning(msg, Fnames{i}, d1, d2);
		testStatus{1,1} = [ pfx.fail msg ];
	    else
		msg = sprintf('CST.Distances.%s = %.6g == %.6g', ...
			       Fnames{i}, d1, d2);
		testStatus{1,1} = [ pfx.ok msg ];
	    end
	    status = [ status ; testStatus ];
	end

    end % isfield(CST, 'Distances')

end % MF version >= 6.0

% Determine if ZEMAX prescription uses segments with equal sized projections.

yVals = CST.ZMNonSeq(2,4,CST.ZMNonSeq(1,4,:)==0 & CST.ZMNonSeq(2,4,:)>0);
yVals = sort(squeeze(yVals));
yVals = diff(yVals);
ratio = max(yVals) / min(yVals);

msg = 'ZEMAX segments have equal projection area';

if ratio < 1.001
    testStatus{1,1} = [ pfx.yes msg ];
else
    testStatus{1,1} = [ pfx.no msg ];
end

status = [ status ; testStatus ];


% Check the vlCsRot* and vlCsEuler routines.

mismatches = zeros(OC.NumSeqs, 1);


% Check the ZEMAX transformations.

for k=1:OC.NumSeqs

    T = CST.ZMSeq(:,:,k);
    C = CST.ZCSeq(:,:,k);

    surfName = OC.SeqSurfName{k};

    if chkZemax 

        % Check the ZEMAX transform against the IM transform.

        Z = vlMzGetGlobalMatrix(k);
        dT = T - Z;
        err = norm(dT) / max(norm(T), norm(Z));

        msg = sprintf('%s: ZMSeq versus ZEMAX (err=%.1g)', surfName, err);

        if err < zmxSeqThreshold
            testStatus{1,1} = [ pfx.ok msg ];
        else
            errMsg = '%s coordinate trans. mismatch CST.ZMSeq(%d) != ZEMAX(%d)';
            warning(errMsg, surfName, k, k)
            mismatches(k) = 1;
            testStatus{1,1} = [ pfx.fail msg ];
        end

        status = [ status ; testStatus ];

    end


    % Compare the CST.ZMSeq transforms with the CST.ZCSeq transforms.
    % Only compare the "rotational" portion of the transform since the
    % support structure is located at a different position than the optics.

    dT = T(:,1:3) - C(:,1:3);
    err = norm(dT);



    % Allow a larger discrepancy between the structural and optical transforms
    % since the structural transforms can only be expected to be "close"
    % in orientation (with an considerable difference in translation values).


    if err < 1.0e-6  % adjust the threshold if required
        msg = sprintf('%s: ZMSeq versus ZCSeq orientation', surfName);
        testStatus{1,1} = [ pfx.ok msg ];
    else
        errMsg = '%s coord. trans. mismatch CST.ZMSeq(%d) != CST.ZCSeq(%d)';
        warning(errMsg, surfName, k, k)
        mismatches(k) = 1;

        msg = sprintf('%s: ZMSeq not aligned with ZCSeq orientation', surfName);
        testStatus{1,1} = [ pfx.warn msg ];
    end

    status = [ status ; testStatus ];


    % Sanity check on CST.CMSeq.

    dT = vlCsMult(C, CST.CMSeq(:,:,k)) - T;
    err = norm(dT) / norm(T);
    msg = sprintf('%s: ZMSeq = ZCSeq * CMSeq check (err=%.2g)', surfName, err);

    if err < 1.0e-15  % adjust the threshold if required
        testStatus{1,1} = [ pfx.ok msg ];
    else
        errMsg = '%s : CST.ZMSeq(%d) != CST.ZCSeq(%d) * CST.CMSeq(%d)';
        warning(errMsg, surfName, k, k, k)
        mismatches(k) = 1;
        testStatus{1,1} = [ pfx.fail msg ];
    end

    status = [ status ; testStatus ];


    % Test alignment of C and M using CMSeq.  This test is similar to the
    % test of CST.ZMSeq versus CST.ZCSeq performed above.

    dT = CST.CMSeq(:,1:3,k) - eye(3);
    err = norm(dT);
    msg = sprintf('%s: Alignment of C%d and M%d', surfName, k+1, k+1);

    if err < 1.0e-6  % adjust the threshold if required
        testStatus{1,1} = [ pfx.ok msg ];
    elseif CST.CMSeq(3,3,k) > 0.999
        mismatches(k) = 1;

        msg = sprintf('%s : Z axes of C%d and M%d are aligned',...
                      surfName, k+1, k+1);
        testStatus{1,1} = [ pfx.ok msg ];
        status = [ status ; testStatus ];

        msg = sprintf('%s : X and Y axes of C%d and M%d are not aligned',...
                      surfName, k+1, k+1);
        warning(msg);
        testStatus{1,1} = [ pfx.info msg ];
    elseif CST.CMSeq(3,3,k) < -0.999
        msg = sprintf('%s : Z axes of C%d and M%d are aligned (opp. dir.)',...
                      surfName, k+1, k+1);
        warning(msg);
        mismatches(k) = 1;
        testStatus{1,1} = [ pfx.info msg ];
    else
        msg = sprintf('%s : Z axes of C%d and M%d are not parallel', ...
                      surfName, k+1, k+1);
        warning(msg);
        mismatches(k) = 1;
        testStatus{1,1} = [ pfx.fail msg ];
    end

    status = [ status ; testStatus ];

    
    % Report offset between C and M using CMSeq.

    dT = abs(CST.CMSeq(3,4,k));
    msg = sprintf('%s: Separation of C%d and M%d = %f', surfName, k, k, dT);

    testStatus{1,1} = [ pfx.info msg ];
    status = [ status ; testStatus ];


    % Print the transformations if there was a discrepancy discovered.

    if mismatches(k)
        disp([ 'CST.ZMSeq   ('  surfName  ')' ] );
        disp(T);

        if chkZemax
            disp([ 'ZEMAX   ('  surfName  ')' ] );
            disp(Z);
        end

        disp([ 'CST.ZCSeq   ('  surfName  ')' ] );
        disp(CST.ZCSeq(:,:,k));
    end
end



% 

if chkZemax 
    Tm1 = vlMzGetGlobalMatrix(0);
end

ConSeg = vlUtSpiralSegtoConSeg((1:OC.NumNonSeqs)');

zemaxErr = zeros(OC.NumNonSeqs,1);
zemaxMaxErr = 0;
cmNorm = zeros(OC.NumNonSeqs,1);
cmErr = [ ];
cmMaxErr = 0;
cmDist = zeros(OC.NumNonSeqs,1);

for k=1:OC.NumNonSeqs


    T = CST.ZMNonSeq(:,:,k);
    C = CST.ZCNonSeq(:,:,k);

    cmNorm(k) = norm(T(:,1:3) - C(:,1:3));

   
    con = ConSeg(k, 1);
    seg = ConSeg(k, 2);

    if chkZemax 

        if all(OC.SegMiss ~= seg)

            Tseg = vlMzGetNSCMatrix(con, seg);
            Tseg = vlCsMult(Tm1, Tseg);

            zemaxErr(k) = norm(T-Tseg) / max(norm(T), norm(Tseg));

            if zemaxErr(k) > zemaxMaxErr
                zemaxMaxErr = zemaxErr(k);
            end;
        end
    end


    % Check the calculation of CST.CMNonSeq

    dT = vlCsMult(C, CST.CMNonSeq(:,:,k)) - T;

    err = norm(dT) / norm(T);
    if err > 1.0e-15
        cmErr = [ cmErr, k ];
    end
    cmMaxErr = max(cmMaxErr, err);

 
    % Store offset between C and M using CMNonSeq, filling in missing segments
    % with -OC.PrimOffset.

    if any(OC.SegMiss == seg)
        % Changed to -OC.PrimOffset (from +OC.PrimOffset) for MF 7.0.
        % For versions prior to 7.0, the cmDist values get converted to
        % absolute values so the change should not affect early versions.

        cmDist(k) = -OC.PrimOffset;
    else
        cmDist(k) = CST.CMNonSeq(3,4,k);
    end
end


% Checks on CMNonSeq transforms

testStatus = cell(1,1);

msg = sprintf('Calculation of CST.CMNonSeq transforms (maxErr=%.2g)', ...
              cmMaxErr);

if isempty(cmErr)
    testStatus{1,1} = [ pfx.ok msg ];
else
    testStatus{1,1} = [ pfx.fail msg ];
    status = [ status ; testStatus ];

    nInvalid = numel(cmErr);

    msg = sprintf('%d invalid CST.CMNonSeq transforms.', nInvalid);
    warning(msg);
    testStatus{1,1} = [ pfx.blank msg ];
    status = [ status ; testStatus ];

    msg = sprintf('First invalid CST.CMNonSeq transform at segment %d', ...
                  cmErr(1));
    warning(msg);
    testStatus{1,1} = [ pfx.blank msg ];
end

status = [ status ; testStatus ];
    

% Report on segment Z axes of C and M which are parallel.
% Each of the missing segments have perfectly aligned coordinate systems
% but we exclude them from the count.

zThreshold = 0.999;

testStatus{1,1} = [ pfx.info ...
    'M1 segment Z axis threshold = ' num2str(zThreshold) ];
status = [ status ; testStatus ];

zNum = sum(CST.CMNonSeq(3,3,:) > zThreshold);

zNum = zNum - 6 * numel(OC.SegMiss); % Subtract missing segments

msg = ' Z axes of C and M segments are aligned.';
testStatus{1,1} = [ pfx.info num2str(zNum) msg ];
status = [ status ; testStatus ];


% Report on segment Z axes of C and M which are anti-parallel.

zNum = sum(CST.CMNonSeq(3,3,:) < -zThreshold);
msg = ' Z axes of C and M segments are aligned (opposite direction).';
testStatus{1,1} = [ pfx.info num2str(zNum) msg ];
status = [ status ; testStatus ];

% Report on segment Z axes of C and M which are not parallel.

npSegs = find(abs(CST.CMNonSeq(3,3,:)) < zThreshold);
zNum = numel(npSegs);

msg = ' Z axes of C and M segments are not parallel.';

if zNum > 0 
    testStatus{1,1} = [ pfx.fail num2str(zNum) msg ];
    status = [ status ; testStatus ];
    testStatus{1,1} = [ pfx.blank 'Segments: ' num2str(npSegs') ];

    status = [ status ; testStatus ];
    testStatus{1,1} = [ pfx.blank 'Minimum Z axes dot product = ' ...
       num2str(min(abs(CST.CMNonSeq(3,3,:)))) ];
else
    testStatus{1,1} = [ pfx.info num2str(zNum) msg ];
end
status = [ status ; testStatus ];



% Results for ZMNonSeq verses ZCNonSeq

cmThreshold = 0.0436;  % Allow 2.5 degree rotation

testStatus = cell(1,1);

if max(cmNorm) <= cmThreshold
    msg = 'Orientation of CST.ZMNonSeq vs. CST.ZCNonSeq';
    testStatus{1,1} = [ pfx.ok msg ];
    status = [ status ; testStatus ];
else
    msg = sprintf('CST.ZMNonSeq & CST.ZCNonSeq are not aligned on %d segments.', ...
                  nnz(cmNorm > cmThreshold));
    testStatus{1,1} = [ pfx.warn msg ];
    status = [ status ; testStatus ];

    errCount = 0;

    for k=1:OC.NumNonSeqs
        if (cmNorm(k) > cmThreshold)
            msg = sprintf('Segment %d: norm(ZMNonSeq - ZCNonSeq) = %g', ...
                          k, cmNorm(k));
            testStatus{1,1} = [ pfx.blank msg ];
            status = [ status ; testStatus ];
            errCount = errCount + 1;
            if errCount >= 20
                break;
            end
        end
    end
end


% Results for M-C separation

if vlSciMeritGetVersion >= 7.0

    % Allow C-M distances to differ from exact value by up to 100 microns.
    cmDistThreshold = 0.000100   % Threshold is 100 microns

    % The cmDist values should be equal to -OC.PrimOffset so we calculate 
    % the error by adding OC.PrimOffset.
    cmDistError = cmDist + OC.PrimOffset;

    % Find the number of C-M distances that are too "small".
    % Since C-M distances are negative, this looks for a positive error term.
    % (corresponding to an absolute C-M distance which is too small).

    numSmall = nnz(cmDistError >= cmDistThreshold);

    % Find the number of C-M distances that are too "large".
    % Since C-M distances are negative, this looks for a negative error term
    % (since the absolute C-M distance is too large).

    numLarge = nnz(cmDistError <= -cmDistThreshold);
else

    cmDist = abs(cmDist);
    numSmall = nnz(abs(cmDist) < 0.5 * OC.PrimOffset);
    numLarge = nnz(abs(cmDist) > 1.00 * OC.PrimOffset);

end

testStatus = cell(1,1);
msg = 'C1-M1 distances for segments'; 

if numSmall + numLarge == 0
    testStatus{1,1} = [ pfx.ok msg ];
    status = [ status ; testStatus ];
else
    testStatus{1,1} = [ pfx.fail msg ];
    status = [ status ; testStatus ];
 
    if numSmall > 0

        % Need to use abs(cmDist) to correct for values being negative.
        % The cmDist values are not converted to absolute values for
        % Merit Function 7.0 and greater.

        msg = sprintf('%d C-M distances are too small (min=%.4f)', ...
                      numSmall, min(abs(cmDist)));
        testStatus{1,1} = [ pfx.blank msg ];
        status = [ status ; testStatus ];
    end

    if numLarge > 0

        % Need to use abs(cmDist) to correct for values being negative.
        % The cmDist values are not converted to absolute values for
        % Merit Function 7.0 and greater.

        msg = sprintf('%d C-M distances are too large (max=%.4f)', ...
                  numLarge, max(abs(cmDist)));
        testStatus{1,1} = [ pfx.blank msg ];
        status = [ status ; testStatus ];

    end

    % For Merit Function 7.0 and greater we print the largest difference
    % between cmDist and -OC.PrimOffset.

    if vlSciMeritGetVersion >= 7.0
        msg = sprintf('Largest C-M discrepancy = %.1e)', ...
                       max(abs(cmDistError)));
        testStatus{1,1} = [ pfx.blank msg ];
        status = [ status ; testStatus ];
    end

end



if chkZemax
    testStatus = cell(1,1);
    msg = sprintf('CST.ZMNonSeq = ZEMAX (maxErr=%.2g)', zemaxMaxErr);

    if max(zemaxErr) <= zmxNonSeqThreshold
        testStatus{1,1} = [ pfx.ok msg ];
        status = [ status ; testStatus ];
    else
        testStatus{1,1} = [ pfx.fail msg ];
        status = [ status ; testStatus ];

        msg = sprintf('CST.ZMNonSeq <> ZEMAX on %d segments.', ...
                      nnz(zemaxErr > zmxNonSeqThreshold));
        testStatus{1,1} = [ pfx.info msg ];
        status = [ status ; testStatus ];

        errCount = 0;
        for k=1:OC.NumNonSeqs
            if (zemaxErr(k) > zmxNonSeqThreshold)
                msg = sprintf('Segment %d: norm(CST.ZMNonSeq - ZEMAX) = %g', ...
                              k, zemaxErr(k));
                testStatus{1,1} = [ pfx.blank msg ];
                status = [ status ; testStatus ];
                errCount = errCount + 1;
                if errCount >= 20
                    break;
                end
            end
        end
    end
end

% Check the ANSYS <--> IM transformation

id = vlSdGetInterColumn('id');
pSegType = vlSdGetIdentNum('PrimarySeg');

xPosId = vlSdGetInterColumn('xPos');
yPosId = vlSdGetInterColumn('yPos');
zPosId = vlSdGetInterColumn('zPos');

ss = [ 'S' num2str(vlSdGetSupportNum(pSegType, 'exists')) ];
existsId = vlSdGetInterColumn(ss);

rowSel = (SD.FEAInterface(:,id) == pSegType);
% nnz(rowSel)

xyzSeg = SD.FEAInterface(rowSel, [ xPosId, yPosId, zPosId, existsId ]);

% size(xyzSeg)

pTriadType = vlSdGetIdentNum('PrimaryTriad');

rowSel = (SD.FEAInterface(:,id) == pTriadType);

ss = [ 'S' num2str(vlSdGetSupportNum(pTriadType, 'segSupport')) ];

segId = vlSdGetInterColumn(ss);

ss = [ 'S' num2str(vlSdGetSupportNum(pTriadType, 'triadNum')) ];

triadId = vlSdGetInterColumn(ss);

xyzTriad = SD.FEAInterface(rowSel, [ xPosId, yPosId, zPosId, segId, triadId ]);


% Find maximum distance between any two segments

segDist = zeros(OC.NumNonSeqs, 1);
mirDist = zeros(OC.NumNonSeqs, 1);
triadDist = zeros(OC.NumNonSeqs, 1);
outerRing = 1 + 6 * sum(1:OC.NumRings-1);

for k=outerRing:OC.NumNonSeqs-3*OC.NumRings;
    k2 = k + 3 * OC.NumRings;
    if xyzSeg(k,4) && xyzSeg(k2,4)
        P1 = xyzSeg(k,1:3)';
        P2 = xyzSeg(k2,1:3)';
        if abs(P1(3) - P2(3)) > 1.0e-6
            warning('...');
        end
        segDist(k) = norm(P1-P2);

        P1 = CST.ZMNonSeq(:,4,k)';
        P2 = CST.ZMNonSeq(:,4,k2)';
        if abs(P1(3) - P2(3)) > 1.0e-6
            warning('...');
        end
        mirDist(k) = norm(P1-P2);

        P1 = CST.ZCNonSeq(:,4,k)';
        P2 = CST.ZCNonSeq(:,4,k2)';
        if abs(P1(3) - P2(3)) > 1.0e-6
            warning('...');
        end
        triadDist(k) = norm(P1-P2);
    end
end


% M1 mirror to mirror distances using CST.ZMNonSeq

msg = sprintf('Maximum distance between M1 mirror centers: %.3f', ...
              max(mirDist));
disp(msg);

testStatus{1,1} = [ pfx.info msg ];
status = [ status ; testStatus ];


% Segment to segment distances using Primary Segment locations from FEA

msg = sprintf('Maximum distance between M1 segment centers: %.3f', ...
              max(segDist));
disp(msg);

testStatus{1,1} = [ pfx.info msg ];
status = [ status ; testStatus ];



% M1 triad to triad distances using CST.ZCNonSeq

msg = sprintf('Maximum distance between M1 triad centers: %.3f', ...
               max(triadDist));
disp(msg);

testStatus{1,1} = [ pfx.info msg ];
status = [ status ; testStatus ];


% Compute the distance between the center of each segment and ...

% delta(i) = norm(vlCsPMult(CST.ZG, xyzSeg_m(k, 1:3)') - vlCsPMult(CST.ZMNonSeq(:,:,k), [ 0 0 0 ]'));

if doGraphs > 0
    st = vlCsValPlots(model, xyzSeg, xyzTriad, doGraphs-1);
end


% end of vlCsCheckTrans function

function status = checkTransformSize(transName, transSize, expSize)

status = [];

pfx = vlUtValPrefixes;

testStatus = cell(1,1);

if numel(transSize) < 3
    msg = sprintf('Size of %s transform (%d, %d)', ...
                  transName, transSize(1), transSize(2));
else
    msg = sprintf('Size of %s transforms (%d, %d, %d)', ...
                  transName, transSize(1), transSize(2), transSize(3));
end

if all(transSize == expSize)
    testStatus{1,1} = [ pfx.ok msg ];
else
    testStatus{1,1} = [ pfx.fail msg ];
    status = [ status ; testStatus ];

    if numel(expSize) < 3
        msg = sprintf('Expected transform size is (%d, %d)', ...
                      expSize(1), expSize(2));
    else
        msg = sprintf('Expected transform size is (%d, %d, %d)', ...
                      expSize(1), expSize(2), expSize(3));
    end
    testStatus{1,1} = [ pfx.blank msg ];

end

status = [ status ; testStatus ];

% end of vlCsCheckTrans.m

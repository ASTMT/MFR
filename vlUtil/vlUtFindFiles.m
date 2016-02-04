% File: <vlUtFindFiles.m>
%
% Syntax: filelist = vlUtilFindFiles(searchdir,matchpattern, recurseflag,
%                       filelist)
%
% Description:
%       Finds all files fitting [matchpattern] in the path [searchdir].
%       Also searches subdirectories if recurseflag is TRUE.  Filelist is a
%       cell array of file names. On input filelist is usually set to {},
%       and is included as an argument for recursion within the function.
%
% Input Parameters:
%       searchdir       - (String) Directory path to search for files
%       matchpattern    - (String) File pattern to match
%       recurseflag     - (binary) If TRUE searches subdirectories to
%                           searchdir
%       filelist        - (Cell Array) List of files matching search criteria. Usually
%                           set as {} on call to function, but used by
%                           function in recursion of subdirs.
%
% Output Parameters:
%       filelist   - (Cell Array) List of files with path matching search criteria.

%
% Required Global Data Structures:
%       None.
%
% Required Data Files:
%       None.
%       

%
% Extended Documentation (Won't be shown in Matlab help command)
%

%
% Revision History
%
% $Id: vlUtFindFiles.m,v 1.1 2006/03/27 17:42:45 roberts Exp $
%
% INDENT-OFF*
% $Log: vlUtFindFiles.m,v $
% Revision 1.1  2006/03/27 17:42:45  roberts
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

function filelist = vlUtilFindFiles(searchdir,matchpattern, recurseflag, filelist)

if nargin ~= 4
    error('Syntax is vlUtilFindFiles(searchdir,matchpattern, recurseflag, {})');
end

% remove any file separator characters from the end of the path string
if searchdir(end) == '\' || searchdir(end) == '/'
    searchdir(end) = [];
end

% make sure were using the correct file separator for the platform
searchdir = strrep(searchdir, '\', filesep);
searchdir = strrep(searchdir, '/', filesep);

pattern = [searchdir filesep matchpattern];
%fprintf('Searching in %s for %s\n',searchdir,pattern);
% add file matches to the list
d = dir(pattern);
for ii = 1:length(d)
    filelist{end+1} = [searchdir filesep d(ii).name];
end

% if recursing look for directories
if recurseflag
    d = dir(searchdir);
    % first two directories will be '.' and '..', so start at 3
    for ii = 3:length(d)
        if d(ii).isdir
            filelist = vlUtFindFiles([searchdir filesep d(ii).name],matchpattern, recurseflag, filelist);
        end
    end
end



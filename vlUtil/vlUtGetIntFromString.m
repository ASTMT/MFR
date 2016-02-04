function [ num ] = vlUtGetIntFromString( str )
% Returns as an integer the first string of numbers found in str, or [] if
% no integer is found

    num = [];

    ii = 1;
    while ~local_isint(str(ii)) 
        if ii == length(str)
            return
        end
        ii = ii+1;
    end

    idxstart = ii;

    while local_isint(str(ii))
        if ii == length(str)
            num = str(idxstart:ii);
            return
        end
        ii = ii+1;
    end
    num = str2num(str(idxstart:ii-1));
end

function isint = local_isint(ch)
    if (ch  >= '0' && ch <= '9')
        isint = 1;
    else
        isint = 0;
    end
end


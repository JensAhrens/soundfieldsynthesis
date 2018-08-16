function [ colors ] = get_colors( N, color_map )
%function [ colors ] = get_colors( N, color_map )
% N                   : desired number of elements
% color_map (optional): string-ID of desired colormap (default: jet)

if ( nargin < 2 )
    colors = jet( N );
else
    function_handle = str2func( color_map );
    colors = function_handle( N );    
end

end


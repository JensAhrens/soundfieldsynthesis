function [ legend_array ] = add_string_to_legend( legend_array, string )
% CREATE_LEGEND_STRING() needs to be called first.
% Finally, call "LEGEND( legend_array )".

legend_array = [ legend_array; string ];

end


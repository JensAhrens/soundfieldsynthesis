function [] = all_labels_latex()
%ALL_LABELS_LATEX Switches all strings in figure to latex
%   Detailed explanation goes here

set( findall( gcf, 'Interpreter', 'tex' ), 'Interpreter', 'Latex' );

end


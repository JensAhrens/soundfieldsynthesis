function [] = picture_size( h, font, font_size, pixels )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ( isempty( h ) )
    h = gca;
end

set( h, 'FontName', font );
set( h, 'FontSize', font_size );
set( gcf, 'Color' , [ 1, 1, 1 ] );  % white background
    
% make sure that the saved image has specified size   
set( gcf, 'PaperUnits', 'inches', 'Position', [ 100 100 pixels pixels] ); 

end


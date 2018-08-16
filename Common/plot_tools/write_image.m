function [ ] = write_image( h, file_name, pixels )
%E.g.: write_image( gcf, file_name, pixels )

if ( nargin < 3 )
    pixels = 200;
end

if ( length( pixels ) == 1 )
    pixels = [ pixels, pixels ];
end

dpi = 256;

% set size
set( h, 'PaperUnits', 'inches', 'Position', [ 100 100 pixels(1) pixels(2) ] );

% white background
set( gcf, 'Color', [ 1, 1, 1 ] );

box on;

% get data and save it
image = getframe( h );
imwrite( image.cdata, file_name, 'png', 'ResolutionUnit', 'meter', 'XResolution', dpi, 'YResolution', dpi );


end


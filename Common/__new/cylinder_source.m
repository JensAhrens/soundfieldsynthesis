function [ S ] = cylinder_source( x, y, x_s, y_s, k )
%function [ S ] = cylinder_source( x, y, x_s, y_s, k )
%
% calculates the sound field at the grid points given by x and y
% of a monochromatic cylinder source located at [ x_s y_s ];
% k = 2*pi*f/c
%

r = sqrt( ( x - x_s ).^2 + ( y - y_s ).^2 );

S = 1i / 4 .* besselh( 0, 2, k .* r );
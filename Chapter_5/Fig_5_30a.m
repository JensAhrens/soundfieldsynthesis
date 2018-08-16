%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This work is supplementary material for the book                        %
%                                                                         %
% Jens Ahrens, Analytic Methods of Sound Field Synthesis, Springer-Verlag %
% Berlin Heidelberg, 2012, http://dx.doi.org/10.1007/978-3-642-25743-8    %
%                                                                         %
% It has been downloaded from http://soundfieldsynthesis.org and is       %
% licensed under a Creative Commons Attribution-NonCommercial-ShareAlike  %
% 3.0 Unported License. Please cite the book appropriately if you use     %
% these materials in your own work.                                       %
%                                                                         %
% (c) 2012 by Jens Ahrens                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% location of point source
x_s = -1;
y_s =  0;

f = 1000;
c = 343;
k = 2*pi*f/c;

% create spatial grid
resolution = 200;
X          = linspace( -2, 2, resolution );
Y          = linspace( -2, 2, resolution );
[ x, y ]   = meshgrid( X, Y );

r = sqrt( ( x - x_s ).^2 + ( y - y_s ).^2 );

% Eq. (1.4)
S = 1/(4*pi) .* exp( -1i .* k .* r ) ./ r;

% normalize to figure center
S = S ./ max( abs( S( end/2, end/2 ) ) );

figure;

imagesc( X, Y, real( S ), [ -2 2 ] );

colormap gray;
turn_imagesc;
axis square;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;

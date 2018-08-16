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

% both Fig. 5.11(a) and (b) will be created 

y_s = -1;
x_s =  0;
f   = 1000; % Hz
c   = 343;  % m/s
k   = (2*pi*f)/c;

% create spatial grid
X        = linspace( -2, 2, 200 );
Y        = linspace( -1, 3, 200 );
[ x, y ] = meshgrid( X, Y );
z        = 0;

% Sec. 3.6.3: Desired field is perfectly synthesized thought mirrored
% at the x-z-plane
r = sqrt( ( x - x_s ).^2 + ( abs( y ) - y_s ).^2 + z.^2 );
S = exp( 1i .* k .* r ) ./ r; % Eq. (1.4)

% normalize
S = S ./ abs( S( end/2, end/2 ) );

% Fig. 5.11(a)
figure;
imagesc( X, Y, real( S ), [ -2 2 ] );
turn_imagesc;
colormap gray;
axis square;

hold on;
% plot secondary source distribution
plot( [ -2 2 ],  [ 0 0 ], 'k', 'Linewidth', 2 )
hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;

% Fig. 5.11(b)
figure;
imagesc( X, Y, 20*log10 ( abs( S ) ), [ -10 10 ] );
turn_imagesc;
colormap gray;
revert_colormap;
colorbar;
axis square;

hold on;
% plot secondary source distribution
plot( [ -2 2 ],  [ 0 0 ], 'k', 'Linewidth', 2 )
hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;

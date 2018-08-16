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

% creates both Fig. 4.42(a) and (b)

clear;

theta_pw = pi/4;
phi_pw   = pi/2;
dx       = 0.2;
y_ref    = 1;
f        = 1300; 
c        = 343;
k        = (2.*pi.*f)./c;

k_pw_x = k .* cos( theta_pw ) .* sin( phi_pw );
k_pw_y = k .* sin( theta_pw ) .* sin( phi_pw );

% create spatial grid
X        = linspace( -2, 4, 500 );
Y        = linspace( -1, 5, 500 );
[ x, y ] = meshgrid( X, Y );

% initialize S
S = zeros( size( x ) );

x_0 = -1 : dx : 1;
y_0 = 0;

% Eq. (3.79)
D = ( 4*1i .* exp( -1i * k_pw_y * y_ref ) ) ./ ( besselh( 0, 2, k_pw_y * y_ref ) ) .* exp( -1i .* k_pw_x .* x_0 );

for index = 1 : length( x_0 )
    S = S + D( index ) .* point_source( x, y, x_0( index ), y_0, 0, k );
end

% normalize
S = S ./ abs( S( round( end/3 ), end/2 ) );

% Fig. 4.42(a)
figure;
imagesc( X, Y, real( S ), [ -2 2 ] );
colormap gray;
turn_imagesc;
axis square

hold on;
% plot secondary source distribution
plot( [ x_0(1) x_0(end) ], [ y_0 y_0 ], 'k', 'LineWidth', 2 );
hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;

% Fig. 4.42(b)

figure;
imagesc( X, Y, 20*log10( abs( S ) ), [ -40 10 ] );
colormap gray;
revert_colormap;
colorbar;
turn_imagesc;
axis square

hold on;
% plot secondary source distribution
plot( [ x_0(1) x_0(end) ], [ y_0 y_0 ], 'k', 'LineWidth', 2 );
hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;

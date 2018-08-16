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

subfigure = 'a'; % for Fig. 4.39(a)
%subfigure = 'b'; % for Fig. 4.39(b)

y = 1;

% frequency axis
f     = linspace( 0, 3000, 200 );
omega = 2.*pi.*f;
c     = 343;

% interval along the x axis that will be sampled
spatial_interval = [ -10 10 ];

% number of spatial sampling points, this makes delta_x = 0.005 m
N = 4001;

delta_x = ( spatial_interval( 2 ) - spatial_interval( 1 ) ) / N;

% this the sampling wavenumber (similar to a sampling frequency f_s)
k_x_s = (2*pi) / delta_x;

% create wavenumber axis; only positive half
k_x    = linspace( 0, k_x_s/2, (N+1)/2 );
k_x(1) = 5*eps; % to avoid numerical instabilities

% now make the negtive part of the wavenumber axis
k_x = [ fliplr( -k_x( 2 : end ) ), k_x ];

% create 2D grid
[ k_x_m omega_m ] = meshgrid( k_x, omega );

% initialize G_kx
G_kx = zeros( size( omega_m ) );

% Eq. (C.10), first case
G_kx( abs( k_x_m ) <= omega_m./c ) = -1i/4 * ...
    besselh( 0, 2, sqrt( ( omega_m( abs( k_x_m ) <= omega_m./c ) ./ c ).^2 - k_x_m( abs( k_x_m ) <= omega_m./c ).^2 ) .* y );

% Eq. (C.10), second case
G_kx( abs( k_x_m ) > omega_m./c ) = 1/(2*pi) * ...
    besselk( 0, sqrt( k_x_m( abs( k_x_m ) > omega_m./c ).^2 - ( omega_m( abs( k_x_m ) > omega_m./c ) ./c ).^2 ) .* y );

if ( subfigure == 'b' )
    % Eq. (4.55)
    G_kx( abs(k_x_m) > 15 ) = 5*eps;
end

figure;
imagesc( k_x, f, 20*log10( abs( G_kx ) ), [ -80 20 ] )

colormap gray; 
revert_colormap;
colorbar;
turn_imagesc;
axis square;

xlim( [ -55 55 ] );

xlabel( 'k_x (rad/m)' );
ylabel( 'f (Hz)' );
graph_defaults;

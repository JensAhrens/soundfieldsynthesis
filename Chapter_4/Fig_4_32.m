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

subfigure = 'a'; % for Fig. 4.32(a)
%subfigure = 'b'; % for Fig. 4.32(b)
%subfigure = 'c'; % for Fig. 4.32(c)
%subfigure = 'd'; % for Fig. 4.32(d)
%subfigure = 'e'; % for Fig. 4.32(e)
%subfigure = 'f'; % for Fig. 4.32(f)

if ( subfigure == 'a' || subfigure == 'c' || subfigure == 'e' )
    f = 1000;
else
    f = 1500;
end

y_ref    = 1;
dx       = .2;
Eta_max   = 5; % highest index of considered spectral repetitions 
c        = 343;
omega    = 2*pi*f;
theta_pw = pi/4;
phi_pw   = pi/2;

% determine which spectral components to consider
if ( subfigure == 'a' || subfigure == 'b' )
    Eta = -Eta_max : Eta_max; % all components 
    
elseif ( subfigure == 'c' || subfigure == 'd' )
    Eta = 0; % only desired component

else
    Eta = [ -Eta_max : -1, 1 : Eta_max ]; % discard desired component

end

% Eq. (A.3)
k_pw_x = omega/c * cos( theta_pw ) * sin( phi_pw );
k_pw_y = omega/c * sin( theta_pw ) * sin( phi_pw );

% create spatial grid
X        = linspace( -2, 2, 300 );
Y        = linspace( -1, 3, 300 );
[ x, y ] = meshgrid( X, Y );
y        = abs( y );

% initialize G_kx and S
G_kx = zeros( size( x ) );
S    = zeros( size( x ) );

% loop over all considered spectral components
for eta = Eta

    % argument of D_tilde in Eq. (4.47)
    k_x = (2*pi)/dx * eta + k_pw_x;
    
    if ( abs(k_x) < omega/c )
        % Eq. (C.10), first case
        G_kx = -1i/4 * besselh( 0, 2, sqrt( (omega/c).^2 - k_x.^2 ) .* y );

    elseif ( abs(k_x) >= omega/c)
        % Eq. (C.10), ssecond case
        G_kx = 1/(2*pi) * besselk( 0, sqrt( k_x.^2 - (omega/c).^2 ) .* y );
           
    end

    % Eq. (4.49)
    S = S + ( 4 * 1i * exp( -i * k_pw_y * y_ref ) ) / ( besselh( 0, 2, k_pw_y .* y_ref ) ) .* ...
            exp( -1i * k_x * x ) .* G_kx;
    
    clear G_kx;
end

figure;
imagesc( X, Y, real( S ), [ -2 2 ] );

hold on;
% plot secondary sources
plot( -2 : dx : 2, 0, 'kx' );
hold off;

turn_imagesc;
colormap gray;
axis square;
xlabel( 'x' );
ylabel( 'y' );
graph_defaults;

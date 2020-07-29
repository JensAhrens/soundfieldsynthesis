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

subfigure = 'a'; % for Fig. 2.7(a)
%subfigure = 'b'; % for Fig. 2.7(b)
%subfigure = 'c'; % for Fig. 2.7(c)
%subfigure = 'd'; % for Fig. 2.7(d)

if ( subfigure == 'a' )
    N = 50;
elseif ( subfigure == 'b' )
    N = 26;
elseif ( subfigure == 'c' )
    N = 13;
elseif ( subfigure == 'd' )
    N = 13;
end

% propagation direction of plane wave
phi_pw   = pi/2; % azimuth
theta_pw = pi/2; % colatitude

f = 1000;
c = 343;
k = 2*pi*f/c;

% create spatial grid
resolution = 300;
X          = linspace( -2, 2, resolution );
Y          = linspace( -2, 2, resolution );
[x, y]     = meshgrid( X, Y );
z          = 0;

% convert to spherical coordinates
r     = sqrt( x.^2 + y.^2 );
alpha = atan2( y, x ); % azimuth
beta  = pi/2;          % colatitude

% initialze S
S = zeros( size( x ) );

% this loop evaluates Eq. (2.38)
for n = 0 : N-1  
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N-1 ) '.' ] );
    
    for m = -n : n
        
        S_breve = 4.*pi .* 1i.^( -n ) .* sphharm( n, -m, theta_pw, phi_pw );
       
        S = S + S_breve .* sphbesselj( n, k.*r ) .* sphharm( n, m, beta, alpha );

    end
end

figure;

if ( subfigure == 'd' )
    imagesc( X, Y, 20*log10( abs( S ) ), [ -20 10 ] );
else
    imagesc( X, Y, real( S ), [ -2 2 ] );
end

colormap gray;

if ( subfigure == 'd' )
    revert_colormap;
    colorbar;
end

turn_imagesc;
axis square;

hold on;

% draw circle with radius r_(N-1), Eq. (2.41)
r_limit  = ( N-1 ) / k;
x_circ   = linspace( -r_limit, r_limit, 200 );

plot( x_circ,  sqrt( r_limit.^2 - x_circ.^2 ), 'k:', 'LineWidth', 2 )
plot( x_circ, -sqrt( r_limit.^2 - x_circ.^2 ), 'k:', 'LineWidth', 2 )

hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;


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

N = 13;

% propagation direction of plane wave
theta_pw = pi/2;
phi_pw   = pi/2;

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
alpha = atan2( y, x );
beta  = pi/2;

% cosine window
window = 0.5 * ( 1 + cos( ( 0 : N ) / N * pi ) );

% initialze S
S = zeros( size( x ) );

% this loop evaluates Eq. (2.38)
for n = 0 : N-1  
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N-1 ) '.' ] );
    
    for m = -n : n
        
        S_breve = 4.*pi .* 1i.^(-n) .* sphharm( n, -m, phi_pw, theta_pw );
        
        S_breve = S_breve .* window( n+1 );

        S = S + S_breve .* sphbesselj( n, k.*r ) .* sphharm( n, m, beta, alpha );
                        
    end
end


figure;

imagesc( X, Y, real( S ), [ -2 2 ] );

colormap gray;
turn_imagesc;
axis square;

hold on;

% draw circle
r_limit  = ( N-1 ) / k;
x_circ   = linspace( -r_limit, r_limit, 200 );

plot( x_circ,  sqrt( r_limit.^2 - x_circ.^2 ), 'k:', 'LineWidth', 2 )
plot( x_circ, -sqrt( r_limit.^2 - x_circ.^2 ), 'k:', 'LineWidth', 2 )

hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;


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

subfigure = 'a';  % for Fig. 4.17(a) 
%subfigure = 'b';  % for Fig. 4.17(b)

if ( subfigure == 'a' )
    M = 27;
else
    M = 79;
end

theta_pw = pi/2;
phi_pw   = pi/2;
L        = 28;
R        = 1.5;
f        = 5000;
c        = 343;
k        = 2*pi*f/c;

% spatial sampling grid
alpha_0 = linspace( 0, 2*pi, 2*L+1 );
alpha_0 = alpha_0( 1:end-1 );
beta_0  = pi/2;

% create spatial grid
X        = linspace( -2, 2, 300 );
Y        = linspace( -2, 2, 300 );
[ x, y ] = meshgrid( X, Y );

% initialize variables
S       = zeros( size( x ) );
G_breve = zeros( 2*M+1, 1 );
S_breve = zeros( 2*M+1, 1 );

for m = -M : M
    % Eq. (2.37a)
    G_breve( m+M+1 ) = -1i .* k .* sphbesselh( abs(m), 2, k.*R ) .* sphharm( abs(m), -m, pi/2, 0 );       
    
    % Eq. (2.38)
    S_breve( m+M+1 ) = 4*pi .* 1i.^( -abs(m) ) .* sphharm( abs(m), -m, phi_pw, theta_pw );
    
end

% loop over secondary sources
for l = 1 : 2*L % loop over all colatitudes

    % initialize D
    D = 0;

    % this loop evaluates Eq. (3.49)
    for m = -M : M    
        D = D + 1 / ( 2*pi*R ) .* S_breve( m+M+1 ) ./ G_breve( m+M+1 ) .* exp( 1i * m* alpha_0( l ) );                
    end

    % position of secondary source
    x_0 = R * cos( alpha_0( l ) ) * sin( beta_0 );
    y_0 = R * sin( alpha_0( l ) ) * sin( beta_0 );
    z_0 = 0;

    S = S + D .* point_source( x, y, x_0, y_0, z_0, k );

end

% normalize 
S = S ./ abs( S( end/2, end/2 ) );

if ( subfigure == 'b' )
    % increase amplitude a bit to account for additional energy due to
    % artifacts
    S = S .* 1.5;
end

figure;

imagesc( X, Y, 20*log10( abs( S ) ), [ -10 10 ] );
turn_imagesc;
colormap gray;
revert_colormap;
colorbar;
axis square;
hold on;

% plot secondary source distribution
plot( R * cos( alpha_0 ), R * sin( alpha_0 ), 'kx' )

if ( subfigure == 'a' )
    % plot r_M region, Eq. (2.41)
    r_limit = M/k;
    x_circ = linspace( -r_limit, r_limit, 100 );
    
    plot( x_circ,  sqrt( r_limit.^2 - x_circ.^2 ), 'k:', 'LineWidth', 2 )
    plot( x_circ, -sqrt( r_limit.^2 - x_circ.^2 ), 'k:', 'LineWidth', 2 )
    
end

hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );

graph_defaults;

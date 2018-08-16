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

subfigure = 'a';  % for Fig. 4.11(a) 
%subfigure = 'b';  % for Fig. 4.11(b)
%subfigure = 'c';  % for Fig. 4.11(c)
%subfigure = 'd';  % for Fig. 4.11(d)
%subfigure = 'e';  % for Fig. 4.11(e)
%subfigure = 'f';  % for Fig. 4.11(f)

if ( subfigure == 'a' || subfigure == 'c' || subfigure == 'e'  )
    N = 28;
else
    N = 79;
end

if ( subfigure == 'a' || subfigure == 'b' )
    f = 1000;
    
elseif ( subfigure == 'c' || subfigure == 'd' )
    f = 2000;
    
elseif ( subfigure == 'e' || subfigure == 'f' )
    f = 5000;
    
end

theta_pw = pi/2;
phi_pw   = pi/2;
L        = 28;
R        = 1.5;
c        = 343;
k        = 2*pi*f/c;

% spatial sampling grid
alpha_0             = linspace( 0, 2*pi, 2*L+1 );
alpha_0             = alpha_0( 1:end-1 );
[ beta_0, weights ] = legpts( L );
beta_0              = acos( beta_0 );

% create spatial grid
X        = linspace( -2, 2, 300 );
Y        = linspace( -2, 2, 300 );
[ x, y ] = meshgrid( X, Y );
z        = 0;

% initialize variables
S       = zeros( size( x ) );
G_breve = zeros( N, 1 );
S_breve = zeros( N, 2*N+1 );

for n = 0 : N-1
    % Eq. (2.37a)
    G_breve( n+1 ) = -1i .* k .* sphbesselh( n, 2, k.*R ) .* sphharm( n, 0, 0, 0 );       
    
    for m = -n : n
        % Eq. (2.38)
        S_breve( n+1, m+N+1 ) = 4*pi .* (-1i).^n .* sphharm( n, -m, phi_pw, theta_pw );
    end
    
end

% loop over secondary sources
for l_1 = 1 : 2*L % loop over all azimuths
    disp( [ 'Calculating azimuth number ' num2str( l_1 ) ' of ' num2str( 2*L ) '.' ] );
    
    for l_2 = 1 : L % loop over all colatitudes
        
        % initialize D
        D = 0;
        
        % this loop evaluates Eq. (3.21)
        for n = 0 : N-1
            for m = -n : n    
                D = D + 1 / ( 2*pi*R^2 ) .* sqrt( ( 2*n+1 ) / (4*pi) ) .* ...
                            S_breve( n+1, m+N+1 ) ./ G_breve( n+1 ) .* ...
                                sphharm( n, m, beta_0( l_2 ), alpha_0( l_1 ) );                
            end
        end
        
        % position of secondary source
        x_0 = R * cos( alpha_0( l_1 ) ) * sin( beta_0( l_2 ) );
        y_0 = R * sin( alpha_0( l_1 ) ) * sin( beta_0( l_2 ) );
        z_0 = R * cos( beta_0( l_2 ) );
        
        S = S + D .* weights( l_2 ) .* point_source( x, y, x_0, y_0, z_0, k );

    end
end

% normalize 
S = S ./ abs( S( end/2, end/2 ) );

figure;

imagesc( X, Y, real( S ), [ -2 2 ] );
turn_imagesc;
colormap gray;
axis square;
hold on;

% plot secondary source distribution
x_circ_0 = linspace( -R, R, 100 );
plot( x_circ_0,  sqrt( R.^2 - x_circ_0.^2 ), 'k--', 'LineWidth', 2 )
plot( x_circ_0, -sqrt( R.^2 - x_circ_0.^2 ), 'k--', 'LineWidth', 2 )

if ( subfigure == 'c' || subfigure == 'e' )
    % plot r_(N-1) region, Eq. (2.41)
    r_limit = ( N-1 )/k;
    x_circ = linspace( -r_limit, r_limit, 100 );
    
    plot( x_circ,  sqrt( r_limit.^2 - x_circ.^2 ), 'k:', 'LineWidth', 2 )
    plot( x_circ, -sqrt( r_limit.^2 - x_circ.^2 ), 'k:', 'LineWidth', 2 )
    
end

hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );

graph_defaults;

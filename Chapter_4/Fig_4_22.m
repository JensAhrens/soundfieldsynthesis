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
    
subfigure = 'a'; % for Fig. 4.22(a)
subfigure = 'b'; % for Fig. 4.22(b)

theta_pw = pi/2;
phi_pw   = pi/2;
L        = 56;
M        = 60;
N_prime  = 15;
R        = 1.5;
c        = 343;
f        = 2000;
k        = 2*pi*f/c;

r_c      = .75;

if ( subfigure == 'a' )
    alpha_c = 0;
elseif ( subfigure == 'b' )
    alpha_c = pi/4;
end

% initialize
S_breve_prime = zeros( N_prime^2, 1 );

for n_prime = 0 : N_prime - 1
    for m_prime = -n_prime : n_prime
        
        % Eq. (2.38)
        S_breve_prime( n_prime^2 + n_prime + m_prime + 1 ) = ...
                                4*pi .* 1i.^( -n_prime ) .* ...
                                            sphharm( n_prime, -m_prime, phi_pw, theta_pw );
                                        
    end
end

% get coefficients (I|I)_(|m|,n')^(m,m'), Eq. (3.55)
II_sect = sectorial_translation_coefficients( M, N_prime, k*r_c, alpha_c, pi/2, 'II' );

% initialize
G_breve = zeros( 2*M+1, 1 );
S_breve = zeros( 2*M+1, 1 );

for m = -M : M
    
    % Eq. (2.37a)
    G_breve( m+M+1 ) = -1i .* k .* sphbesselh( abs(m), 2, k.*R ) .* sphharm( abs(m), -m, pi/2, 0 );       
    
    % this loop evaluates Eq. (3.54)
    for n_prime = 0 : N_prime - 1
        for m_prime = -n_prime : n_prime

            S_breve( m+M+1 ) = S_breve( m+M+1 ) + S_breve_prime( n_prime^2 + n_prime + m_prime + 1 ) * ...
                                                        (-1)^( abs(m) + n_prime ) * ...
                                                            II_sect( n_prime^2 + n_prime + m_prime + 1, m + M + 1 );
            
        end
    end
    
end

% spatial sampling grid
alpha_0 = linspace( 0, 2*pi, L+1 );
alpha_0 = alpha_0( 1 : end-1 );
beta_0  = pi/2;

% create spatial grid
X        = linspace( -2, 2, 300 );
Y        = linspace( -2, 2, 300 );
[ x, y ] = meshgrid( X, Y );

% initialize
S        = zeros( size( x ) );

% loop over secondary sources
for l = 1 : L

    % initialize
    D = 0;

    % this loop evaluates Eq. (3.49)
    for m = -M : M    
        D = D + 1 / ( 2*pi*R ) .* S_breve( m+M+1 ) ./ G_breve( m+M+1 ) .* exp( 1i * m * alpha_0( l ) );                
    end

    % position of secondary source
    x_0 = R * cos( alpha_0( l ) ) * sin( beta_0 );
    y_0 = R * sin( alpha_0( l ) ) * sin( beta_0 );
    z_0 = 0;

    S = S + D .* point_source( x, y, x_0, y_0, z_0, k );

end

figure;

imagesc( X, Y, real( S ), [ -100 100 ] );
turn_imagesc;
colormap gray;
axis square;

hold on;

% plot secondary source distribution
plot( R * cos( alpha_0 ), R * sin( alpha_0 ), 'kx' )

% plot local center
plot( r_c * cos( alpha_c ), r_c * sin( alpha_c ), 'ko', 'Color', [ .99 .99 .99 ] );

hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );

graph_defaults;

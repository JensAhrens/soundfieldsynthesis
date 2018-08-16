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

N = 8;
Q = 18;

r_s     = 3;
alpha_s = pi/4;
beta_s  = pi/2;
 
f = 500;
c = 343;
k = 2*pi*f/c;

% create spatial grid
X        = linspace( -3, 3, 300 );
Y        = linspace( -3, 3, 300 );
[ x, y ] = meshgrid( X, Y );

x_0 = [ linspace( -2, 2, 21 ).'; ...
        2 * ones( 19, 1 ); ...
       linspace( -2, 2, 21 ).'; ...
       -2 * ones( 19, 1 ) ];

y_0 = [ 2 .* ones( 21, 1 ); ...
       linspace( -1.8, 1.8, 19 ).'; ...
       -2 * ones( 21, 1 ); ...
       linspace( -1.8, 1.8, 19 ).' ];
   
alpha_n = [ -pi/2 .* ones( 21, 1 ); ...
             pi   .* ones( 19, 1 ); ...
             pi/2 .* ones( 21, 1 ); ...
             zeros( 19, 1 ) ];
         
% initialize S
S = zeros( size( x ) );

% loop over secondary sources
for l = 1 : length( x_0 )
    
    disp( [ 'Processing secondary source ' num2str( l ) ' of ' num2str( length( x_0 ) ) '.' ] );

    D = 0;

    % this loop evaluates Eq. (5.91)
    for n = 0 : N-1
        for m_prime = 0 : n
 
            if ( abs( -( 2*m_prime - n ) ) <= n )
                
                % Eq. (2.37a)
                S_breve = -1i.* k .* sphbesselh( n, 2, k.*r_s ) .* ...
                                        sphharm( n, -( 2*m_prime - n ), beta_s, alpha_s );
            
                % Eq. (5.91)
                D = D + Psi_projection( n, 2*m_prime - n ) .* 1i^n .* ...
                        S_breve .* Lambda_projection( 2*m_prime - n, alpha_n(l), x_0(l), y_0(l), k, Q );
            end
            
        end
    end

    S = S + D .* point_source( x, y, x_0(l), y_0(l), 0, k );

end

% normalize
S = S ./ abs( S( end/2, end/2 ) );

% rotate sound field by 45 deg to make it look exactly like in Fig. 5.57(c)
S = S .* 1i;

figure;
imagesc( X, Y, real( S ), [ -2 2 ] );
turn_imagesc;
colormap gray;

hold on;

% plot secondary sources
plot( x_0, y_0, 'kx' );

% plot r_(N-1)-boundary
r_limit = (N-1) / k;
x_circ = linspace( -r_limit, r_limit, 200 );

plot( x_circ,  sqrt( r_limit.^2 - x_circ.^2 ), 'k:' )
plot( x_circ, -sqrt( r_limit.^2 - x_circ.^2 ), 'k:' )

hold off;

axis square;
xlabel( 'x (m)' );
ylabel( 'y (m)' );

graph_defaults;

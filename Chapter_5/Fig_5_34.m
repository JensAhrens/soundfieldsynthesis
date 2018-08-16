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

subfigure = 'a'; % for Fig. 5.34(a)
%subfigure = 'b'; % for Fig. 5.34(b)

if ( subfigure == 'a' )
    r_s     = 1;
    alpha_s = pi;
else
    r_s     = 0.5;
    alpha_s = pi/4;
end

R       = 1.5;  % m
f       = 1000; % Hz
c       = 343;  % m/s

k = (2*pi*f)/c;

% from Eq. (5.40)
N = floor( k*r_s );

% create spatial grid
X        = linspace( -2, 2, 200 );
Y        = linspace( -2, 2, 200 );
[ x, y ] = meshgrid( X, Y );

% convert to spherical coordinates
r     = sqrt( x.^2 + y.^2);
alpha = atan2( y, x );
beta  = pi/2;

% initialize S_int and S_ext
S_int = zeros( size( r ) );
S_ext = zeros( size( r ) );

% this loop evaluates Eq. (3.50) and (3.51) with (5.40)
for n = 0 : N
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N ) '.' ] );
       
    for m = -n : n
        
        % Eq. (2.37a)
        S_breve = -1i .* k .* sphbesselh( abs(m), 2, k.*r_s ) .* sphharm( abs(m), -m, pi/2, alpha_s );   
        
        % Eq. (5.40)
        S_breve = S_breve .* 0.5 * ( cos( n / N * pi ) + 1 );
        
        % Eq. (2.37a)
        G_breve_int_nm = -1i .* k .* sphbesselh( n, 2, k.*R ) .* sphharm( n, -m, pi/2, 0 );

        % Eq. (2.37a)
        G_breve_int_mm = -1i .* k .* sphbesselh( abs(m), 2, k.*R ) .* sphharm( abs(m), -m, pi/2, 0 );

        % Eq. (2.37b)
        G_breve_ext    = -1i .* k .* sphbesselj( n, k.*R ) .* sphharm( n, -m, pi/2, 0 );
    
        % Eq. (3.50)
        S_int = S_int + S_breve .* G_breve_int_nm ./ G_breve_int_mm .* ...
                                sphbesselj( n, k.*r ) .* sphharm( n, m, beta, alpha );
        
        % Eq. (3.51)
        S_ext = S_ext + S_breve .* G_breve_ext ./ G_breve_int_mm .* ...
                                sphbesselh( n, 2, k.*r ) .* sphharm( n, m, beta, alpha );      
        
     end
end

% assemble exterior and interior field
S          = zeros( size( x ) ); 
S( r < R ) = S_int( r < R );
S( r > R ) = S_ext( r > R );

% normalize
S = S ./ abs( S( end/2, end/2 ) );

figure;
imagesc( X, Y, real( S ), [ -2 2 ] );
turn_imagesc;
colormap gray;
axis square;

hold on;

% plot secondary source distribution
x_circ = linspace( -R, R, 200 );
plot( x_circ,  sqrt( R.^2 - x_circ.^2 ), 'k', 'Linewidth', 2 )
plot( x_circ, -sqrt( R.^2 - x_circ.^2 ), 'k', 'Linewidth', 2 )

% mark focused source
plot( r_s * cos( alpha_s ), r_s * sin( alpha_s ), 'kx'  )

% plot dotted line
if ( subfigure == 'a' )
    plot( [ -r_s -r_s ], sqrt( R.^2 - r_s.^2 ) .* [ -1 1 ], 'k:' )
else
    plot( [ R * cos( acos(1/3) + pi/4 ), R * sin( acos(1/3) + pi/4 ) ], ...
          [ R * cos( pi/4 - acos(1/3) ), R * sin( pi/4 - acos(1/3) ) ], 'k:' )
end

hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );

graph_defaults;

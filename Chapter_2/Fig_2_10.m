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

subfigure = 'a'; % for Fig. 2.10(a)
%subfigure = 'b'; % for Fig. 2.10(b)

if ( subfigure == 'a' )
    N = 13;
elseif ( subfigure == 'b' )
    N = 26;
end

% location of point source
alpha_s = -pi/2;
beta_s  = pi/2;
r_s     = 1;

f = 1000;
c = 343;
k = 2*pi*f/c;

% create spatial grid
resolution = 300;
X          = linspace( -2, 2, resolution );
Y          = linspace( -2, 2, resolution );
[ x, y ]   = meshgrid( X, Y );
z          = 0;

% convert to spherical coordinates
r     = sqrt( x.^2 + y.^2 );
alpha = atan2( y, x );
beta  = pi/2;

% initialze S
S = zeros( size( x ) );

% this loop evaluates Eq. (2.37a)
for n = 0 : N-1  
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N-1 ) '.' ] );
    
    for m = -n : n
        
        S_breve = -1i .* k .* sphbesselh( n, 2, k*r_s ) .* ...
                                sphharm( n, -m, beta_s, alpha_s );

        S = S + S_breve .* sphbesselj( n, k.*r ) .* sphharm( n, m, beta, alpha );

    end
end

% normalize to figure center
S = S ./ max( abs( S( end/2, end/2) ) );

figure;

imagesc( X, Y, real( S ), [ -2 2 ] );

colormap gray;
turn_imagesc;
axis square;

hold on;

% draw circles
x_circ_1 = linspace( -r_s, r_s, 200 );
r_limit  = ( N-1 ) / k;
x_circ_2 = linspace( -r_limit, r_limit, 200 );

% region of validity of the expansion
plot( x_circ_1,  sqrt( r_s.^2 - x_circ_1.^2 ), 'k--', 'LineWidth', 2)
plot( x_circ_1, -sqrt( r_s.^2 - x_circ_1.^2 ), 'k--', 'LineWidth', 2)
    
if ( subfigure == 'a' )
    % r_(N-1) region, Eq. (2.41)
    plot( x_circ_2,    sqrt( r_limit.^2 - x_circ_2.^2 ), 'k:', 'LineWidth', 2 )
    plot( x_circ_2,   -sqrt( r_limit.^2 - x_circ_2.^2 ), 'k:', 'LineWidth', 2 )
    
    annotation( 'arrow', [ .515 .515 ], [ .17 .27 ] );
    annotation( 'arrow', [ .515 .515 ], [ .77 .87 ] );
    annotation( 'arrow', [ .655 .735 ], [ .73 .85 ] );
    annotation( 'arrow', [ .375 .295 ], [ .73 .85 ] );
end

hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;


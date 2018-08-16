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

subfigure = 'a'; % for Fig. 4.24(a)
%subfigure = 'b'; % for Fig. 4.24(b)

N = 79;
M = 27;
R = 1.5;
c = 343;

if ( subfigure == 'a' )
    f = 1000; 
else
    f = 3000;
end

k = (2*pi*f)/c;

% create spatial grid
X        = linspace( -2, 2, 200 );
Y        = linspace( -2, 2, 200 );
[ x, y ] = meshgrid( X, Y );

r        = sqrt( x.^2 + y.^2 );
alpha    = atan2( y, x );
beta     = pi/2;

% initialize S 
S = zeros( size( r ) );

for n = 0 : N-1
   disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N-1 ) '.' ] );
   
   for m = max( -n, -M ) : min( n, M ) % this results from a combination of Eq. (4.40) and Eq. (2.34)
        
       % Eq. (2.37a)
       S( r<R ) = S( r<R ) - 1i .* k .* sphbesselh( n, 2, k.*R ) .* ...
                    sphharm( n, -m, pi/2, 0 ) .* ...
                        sphbesselj( n, k.*r( r<R ) ) .* ...
                            sphharm( n, m, beta, alpha( r<R ) );

       % Eq. (2.37b) 
       S( r>R ) = S( r>R ) - 1i .* k .* sphbesselh( n, 2, k.*r( r>R ) ) .* ...
                    sphharm( n, -m, pi/2, 0 ) .* ...
                        sphbesselj( n, k.*R ) .* ...
                            sphharm( n, m, beta, alpha( r>R ) );
    
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
% plot position of source
plot( R, 0, 'x', 'Color', [ .1 .1 .1 ] );
hold off;

xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;

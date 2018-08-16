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

subfigure = 'a'; % for Fig. 5.18(a)
%subfigure = 'b'; % for Fig. 5.18(b)

if ( subfigure == 'a' )
    f = 1000;

elseif ( subfigure == 'b' )
    f = 3000;

end

omega = 2*pi*f;
c     = 343;
k     = omega/c;
N     = 31;

L       = 6;
d_x0    = .1;
d_ref   = 1;

alpha_o = pi/3;
beta_o  = pi/2;
alpha_n = pi/2;

x_0 = -L/2 : d_x0 : L/2;
y_0 = 1;

% create spatial grid
X        = linspace( -1, 3, 400 ); 
Y        = linspace(  0, 4, 400 ); 
[ x, y ] = meshgrid( X, Y );

% initialize S
S = zeros( size( x ) );    
    
% loop over secondary sources
for l = 1 : length( x_0 )
    
    disp( [ 'Processing secondary source ' num2str( l ) ' of ' num2str( length( x_0 ) ) '.' ] );
    
    r_0     = sqrt( x_0( l ).^2 + y_0.^2 );
    alpha_0 = atan2( y_0, x_0( l ) );
   
    % initialize D
    D = 0;
    
    % this loop evaluates Eq. (5.18)
    for n = 0 : N-1
        for m = -n : n
            
            % Eq. (2.44)
            S_breve = 1i^( -n ) .* ...
                    ( factorial( N - 1 ) * factorial( N ) ) ./ ...
                        ( factorial( N + n ) * factorial( N - n - 1 ) ) .* ...
                            sphharm( n, -m, beta_o, alpha_o );
            
            % Eq. (5.18)         
            D = D + sqrt( (2*pi*d_ref) / (1i*omega) ) .* ...
                        cos( alpha_n - alpha_0 ) .* exp( -1i .* k .* r_0 ) ./ r_0 .* ...
                            1i^n .* S_breve .* sphharm( n, m, pi/2, alpha_0 );
        end
    end
    
    S = S + D .* point_source( x, y, x_0(l), y_0, 0, k );    
    
end

% normalization
S = S ./ abs( S( end/2, end/2 ) );

figure;
imagesc( X, Y, real( S ), [ -2.5 2.5 ] )
colormap gray;
xlabel( 'x (m)' );
ylabel( 'y (m)' );
turn_imagesc;
axis square;

% plot secondary sources
hold on;
plot( x_0, y_0 .* ones( size( x_0 ) ), 'k.', 'Linewidth', 2 );
hold off;

graph_defaults;

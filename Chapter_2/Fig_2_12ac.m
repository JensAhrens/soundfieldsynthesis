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

subfigure = 'a'; % for Fig. 2.12(a)
%subfigure = 'c'; % for Fig. 2.12(c)

if ( subfigure == 'a' )
    N = 4;
elseif ( subfigure == 'c' )
    N = 21;
end    

% orientation of source
alpha_or = 0;
beta_or  = pi/2;

f = 1000;
c = 343;
k = 2*pi*f/c;

% create spatial grid
X        = linspace( -2, 2, 200 );
Y        = linspace( -2, 2, 200 );
[ x, y ] = meshgrid( X, Y );
z        = 0;

% transfer to spherical coordinates
r     = sqrt( x.^2 + y.^2 );
alpha = atan2( y, x );
beta  = pi/2;

% initliaize S
S = zeros( size( x ) );

% this loop evaluates Eq. (2.32b) with (2.44)
for n = 0 : N-1
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N-1 ) '.' ] );
    
    for m = -n : n
        S_breve = 1i^( -n ) .* ...
            ( factorial( N-1 ) * factorial( N ) ) / ( factorial( N+n ) * factorial( N-n-1 ) ) .* ...
                        sphharm( n, -m, beta_or, alpha_or );
                    
        S = S + S_breve .* sphbesselh( n, 2, k.*r ) .* sphharm( n, m, beta, alpha );
    end
end


% normalize
if ( subfigure == 'a' )
    norm_factor = 20; 
    
elseif ( subfigure == 'c' )
    norm_factor = 10;

end  

S = S .* norm_factor;

figure;
imagesc( X, Y, real( S ), [ -1 1 ] );
turn_imagesc;
colormap gray;
axis square;
xlabel( 'x (m)' );
ylabel( 'y (m)' );
graph_defaults;

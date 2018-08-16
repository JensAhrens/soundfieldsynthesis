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

% Set this number lower if you run out of memory, e.g., to 100.
no_of_frequencies = 300; 

N = 6;

% propagation direction of plane wave
theta_pw = pi/2;

% frequencies to display
f = linspace( 0, 5000, no_of_frequencies ).';
c = 343;
K = 2*pi*f/c;

% create grid
X           = linspace( -2, 2, 200 );
Y           = linspace( -2, 2, 200 );
[ x, y, k ] = meshgrid( X, Y, K );
r           = sqrt( x.^2 + y.^2 );
alpha       = atan2( y, x );

% initialize S
S = zeros( length( X ), length( Y ), length( K ) );

% this loop evaluates Eq. (2.38) for all frequencies in vector f
for n = 0 : N-1
    disp( [ 'Calculating order ' num2str( n ) ' of ' num2str( N-1 ) '.' ] );
    
    for m = -n : n
        S = S + 4 .* pi .* 1i^( -n ) .* sphharm( n, -m, pi/2, theta_pw ) .* ...
                    sphbesselj( n, k.*r ) .* sphharm( n, m, pi/2, alpha );
        
    end
end

% clear some memory and assure access to Matlab function 'alpha'
clear alpha r x y k;

% prepare data for plotting
data = 20*log10( abs( S ) );

clear S;

% clip data
data_min = -30;
data_max =  10;

data( data < data_min ) = data_min;
data( data > data_max ) = data_max;

figure;

% plot data
h = slice( X, Y, f, data, [], [], f );

axis normal;

% create transparency
alpha( 'color' )
set( h, 'EdgeColor', 'none', 'FaceAlpha', 'flat' );

daspect( [ .004 0.004 5 ] )
view( 15, 30 );

caxis( [ data_min data_max ] );
zlim( [ f(1) f(end) ] );

colormap gray;
revert_colormap;
colorbar;

xlabel( '$x$~(m)' , 'FontSize', 13, 'Interpreter', 'latex' );
ylabel( '$y$~(m)' , 'FontSize', 13, 'Interpreter', 'latex' );
zlabel( '$f$~(Hz)', 'FontSize', 13, 'Interpreter', 'latex' );

graph_defaults;



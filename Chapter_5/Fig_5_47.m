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

subfigure = 'a'; % for Fig. 5.47(a)
%subfigure = 'b'; % for Fig. 5.47(b)

if ( subfigure == 'a' )
    v = 120; 
elseif ( subfigure == 'b' )
    v = 240;
end

f     = 500; 
t     = 0;
c     = 343; 
omega = 2*pi*f;
M     = v / c;
dx_0  = 0.1;
d_ref = 1;

% create spatial grid
X        = linspace( -3, 3, 500 );
Y        = linspace( -1, 5, 500 );
[ x, y ] = meshgrid( X, Y );
z        = 0;

% initialize s
s = zeros( size( x ) );

y_0 = 1;
z_0 = 0;

% loop over secondary sources (Eq. (5.69))
for x_0 = -5 : dx_0 : 5

    disp( [ 'Processing secondary source at x_0 = ' num2str( x_0, '%0.1f' ) ' m.' ] ); 
    
    r_0 = sqrt( ( x - x_0 ).^2 + ( y - y_0 ).^2 + ( z - z_0 ).^2  );
    
    % Eq. (5.58) and (5.69)
    Delta = sqrt( ( x_0 - v * ( t - r_0/c ) ).^2  + ( y_0.^2 + z_0.^2 ) .* ( 1 - M^2 ) );

    % Eq. (5.57) and (5.69)
    tau = ( M .* ( x_0 - v * ( t - r_0/c ) ) + Delta ) / ( c * ( 1 - M^2 ) );

    % Eq. (5.66), (5.69), and (3.93)
    d = sqrt( 2*pi*d_ref / ( 1i * omega/c ) ) .* 2 * y_0 ./ Delta .* ( ( 1 - M^2 ) ./ Delta + 1/c * 1i*omega ) .* exp( 1i .* omega .* ( t - r_0/c - tau ) ) ./ Delta;
    s = s + 1 / (4*pi) .* d ./ r_0;

end
    
% normalize
s = s ./ abs( s( end/2, end/2 ) );

figure;
imagesc( X, Y, real( s ), [ -1.5 1.5 ] );
turn_imagesc;
axis square;
colormap gray;
revert_colormap;

hold on;
% plot secondary sources
plot( ( -3 : dx_0 : 3 ), y_0, 'kx' )
% plot position of virtual source
plot( 0, 0, 'k.' )
hold off;

xlabel( 'x (m)' )
ylabel( 'y (m)' )

graph_defaults;

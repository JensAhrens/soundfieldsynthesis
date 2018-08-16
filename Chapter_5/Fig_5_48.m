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

subfigure = 'a'; % for Fig. 5.48(a)
%subfigure = 'b'; % for Fig. 5.48(b)
%subfigure = 'c'; % for Fig. 5.48(c)
%subfigure = 'd'; % for Fig. 5.48(d)

v     = 600; 
f     = 500; 
t     = 0;
c     = 343; 
omega = 2*pi*f;
M     = v / c;
dx_0  = 0.1;
d_ref = 1;

% create spatial grid
X        = linspace( -5, 1, 500 );
Y        = linspace( -1, 5, 500 );
[ x, y ] = meshgrid( X, Y );
z        = 0;

% initialize s
s = zeros( size( x ) );

y_0 = 1;
z_0 = 0;

% loop over secondary sources (Eq. (5.69))
for x_0 = -7 : dx_0 : 0

    disp( [ 'Processing secondary source at x_0 = ' num2str( x_0, '%0.1f' ) ' m.' ] ); 
    
    r_0 = sqrt( ( x - x_0 ).^2 + ( y - y_0 ).^2 + ( z - z_0 ).^2 );
    
    % Eq. (5.58) and (5.69)
    Delta = sqrt( ( x_0 - v * ( t - r_0/c ) ).^2  + ( y_0.^2 + z_0.^2 ) .* ( 1 - M^2 ) );

    % Eq. (5.57) and (5.69)
    tau_1 = ( M .* ( x_0 - v * ( t - r_0/c ) ) + Delta ) / ( c * ( 1 - M^2 ) );

    % Eq. (5.57) and (5.69)
    tau_2 = ( M .* ( x_0 - v * ( t - r_0/c ) ) + Delta ) / ( c * ( 1 - M^2 ) );

    % Eq. (5.66) and (5.69)
    d = sqrt( 2*pi*d_ref / ( 1i * omega/c ) ) * 2 * y_0 ./ Delta .* ( ( 1 - M^2 ) ./ Delta + 1/c * 1i*omega ) .* exp( 1i .* omega .* ( t - r_0/c - tau_1 ) ) ./ Delta;
    
    % Eq. (5.71), (5.69), and (3.93)
    d = d + sqrt( 2*pi*d_ref / ( 1i * omega/c ) ) * 2 * y_0 ./ Delta .* ( ( 1 - M^2 ) ./ Delta - 1/c * 1i*omega ) .* exp( 1i .* omega .* ( t - r_0/c - tau_2 ) ) ./ Delta;
    
    % time instant when Mach cone passes secondary source             
    t_mach = 1/v * ( x_0 + y_0 * sqrt( M^2 - 1 ) );
    
    if ( subfigure == 'a' || subfigure == 'b' )
        t_fade_in = 0; % start of fade-in after Mach cone
    elseif ( subfigure == 'c' )
        t_fade_in = 0.0001; % start of fade-in after Mach cone
    elseif ( subfigure == 'd' )
        t_fade_in = 0.0002; % start of fade-in after Mach cone
    end
        
    if ( subfigure ~= 'a' )
        fade_dur = 0.0001; % duration of fade-in
        
        % fade-in
        d( ( t - r_0/c ) < t_mach + t_fade_in + fade_dur ) = d( ( t - r_0/c ) < t_mach + t_fade_in + fade_dur ) .* ...
                                        sin( ( t - r_0( ( t - r_0/c ) < t_mach + t_fade_in + fade_dur )/c - ( t_mach + t_fade_in ) ) / fade_dur * pi/2 ).^2;
    end
 
    % driving function is zero before Mach cone arrives
    d( ( t - r_0/c ) < t_mach + t_fade_in ) = 0;

    s = s + 1 / (4*pi) .* d ./ r_0;

end
    
figure;
imagesc( X, Y, real( s ), [ -20 20 ] );
turn_imagesc;
axis square;
colormap gray;

hold on;
% plot secondary sources
plot( ( -5 : dx_0 : 1 ), y_0, 'kx' )
hold off;

xlabel( 'x (m)' )
ylabel( 'y (m)' )

graph_defaults;

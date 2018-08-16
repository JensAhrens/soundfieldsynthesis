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

% Fig. 4.19(b), (d), and (f) will be created; we use WFS for convenience

theta_pw = pi/2;
L        = 56;
r_0      = 1.5;
    
no_of_bins = 1024;

fs   = 44100;
c    = 343;

d_alpha_0 = 2 * pi / L;

% time instances to plot
t_to_show    = [ -0.0027, 0, 0.0027 ]; % in sec
taps_to_show = round( t_to_show * fs + 1280 );

maximum_delay = ceil( r_0 / c * fs ); % in samples

d = zeros( maximum_delay + 5 * no_of_bins + 2, L );

% Eq. (3.102)
prefilter = wavread( 'wfs_prefilter_100_1800_44100.wav' );

% loop over secondary sources
for l = 1 : L 
    
    alpha_0 = l * d_alpha_0;
    alpha_n = alpha_0 + pi;
    
    % Eq. (3.89)
    if ( cos( alpha_n - theta_pw ) < 0 ) 
        continue;
    end

    % from Eq. (3.101)
    delay = r_0 / c * ( 1 - cos( alpha_n - theta_pw ) );

    delay = round( delay * fs );
    
    % from Eq. (3.101)
    amplitude = cos( theta_pw - alpha_n );
    
    d( delay + 1 : delay + length( prefilter ), l ) = amplitude .* prefilter;
       
end

% put zeros around to have some headroom
d = [ zeros( no_of_bins, L ); d; zeros( no_of_bins, L ) ];

% normalize
d = d ./ max(abs( d( : ) ) );

% create spatial grid
resolution = 400;
X          = linspace( -2, 2, resolution );
Y          = linspace( -2, 2, resolution );
[ x, y ]   = meshgrid( X, Y );

for tap = taps_to_show

    s          = zeros( size( x ) );
    d_reshaped = zeros( size( x ) );

    % this loop evaluates Eq.(5.69)
    for l = 1 : L

        x_0 = r_0 * cos( l * d_alpha_0 );
        y_0 = r_0 * sin( l * d_alpha_0 );

        r = sqrt( ( x - x_0 ).^2 + ( y - y_0 ).^2 );

        % from Eq. (5.69)
        t = ( tap/fs - r./c ) .* fs + 1; % in samples

        % Interpolate the impulse responses to find the values at instances
        % t, which correspond to the spatial locations that we are 
        % interested in. 
        d_reshaped = reshape( ...
                          interp1( ( 1 : size( d, 1 ) ), d( :, l ), t, 'linear' ), ...
                                                                resolution, resolution );

        % from Eq. (5.69)
        s( find( r > 0 ) ) = s( find( r > 0 ) ) + d_reshaped( find( r > 0 ) ) ./ r( find( r > 0 ) );

    end

    figure;
    % avoid log of 0
    imagesc( X, Y, 20*log10( abs( s ) + 5*eps ), [ -30 0 ] );
    
    hold on;
    
    % plot secondary sources
    plot( r_0 .* cos( ( 0 : L-1 ) .* d_alpha_0 ), r_0 .* sin( ( 0 : L-1 ) .* d_alpha_0 ), 'kx' );
    
    hold off;
    
    axis square;
    turn_imagesc;
    colormap gray;
    revert_colormap;
    colorbar;
    xlabel( 'x (m)' );
    ylabel( 'y (m)' );
    drawnow;

    graph_defaults;

end % tap

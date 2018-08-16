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

subfigure = 'b'; % for Fig. 4.20(b)
%subfigure = 'd'; % for Fig. 4.20(d)

theta_pw = pi/2;
L        = 56;
R        = 1.5;
x        = 1;
y        = 0;

no_of_bins = 1024;

fs     = 44100;
c      = 343;
cutoff = 2000; % for Fig. 4.20(d)

% time instances to plot
t_to_show    = [ -0.004 0.01 ]; % in sec
taps_to_show = round( ( t_to_show(1) * fs : t_to_show(2) * fs ) + 1060 );
    
d_alpha_0 = 2 * pi / L;

maximum_delay = ceil( R / c * fs ); % in samples

d = zeros( maximum_delay + 5 * no_of_bins + 2, L );

prefilter = wavread( 'wfs_prefilter_100_1800_44100.wav' );

% loop over secondary sources
for l = 1 : L 
    
    alpha_0 = l * d_alpha_0;
    alpha_n = alpha_0 + pi;
    
    % Eq. (3.89)
    if ( cos( alpha_n - theta_pw ) < 0 ) 
        continue;
    end

    delay = R / c * ( 1 - cos( alpha_n - theta_pw ) );

    delay = round( delay * fs );
    
    amplitude = cos( theta_pw - alpha_n );
    
    d( delay + 1 : delay + length( prefilter ), l ) = amplitude .* prefilter;
       
end

% put zeros around to have some headroom
d = [ zeros( no_of_bins, L ); d; zeros( no_of_bins, L ) ];

% normalize
d = d ./ max( abs( d( : ) ) );

% initialize s
s = zeros( length( taps_to_show ), 1 );

% this loop evaluates Eq. (5.69)
for l = 1 : L
    
    x_0 = R * cos( (l-1) * d_alpha_0 );
    y_0 = R * sin( (l-1) * d_alpha_0 );
    
    r = sqrt( ( x - x_0 ).^2 + ( y - y_0 ).^2 );
    
    % from Eq. (5.69)
    t = ( taps_to_show./fs - r./c ) .* fs + 1; % in samples

    % Interpolate the impulse responses to find the values at instances
    % t that we are interested in. 
    d_interp = interp1( ( 1 : size( d, 1 ) ), d( :, l ), t, 'linear' ).';
    
    % from Eq. (5.69)
    s = s + d_interp ./ r;
    
end

figure;

if ( subfigure == 'b' )
    % avoid log of 0
    s = 20*log10( abs( s ) + 5*eps );
    
    plot( ( taps_to_show - 1302 ) ./ fs .* 1000, s, 'k' )
    
else
    % apply filtering in FFT domain (ugly but serves its purpose...)
    S   = fft( s );
    S_1 = S;
    S_2 = S;
    
    % cutoff frequency in samples
    cutoff = round( cutoff./fs .* size( s, 1 ) ) + 1; % samples

    % set all portions that we want to suppress to zero
    S_1( 1 : cutoff )           = 0; 
    S_1( end-cutoff+1 : end )   = 0;
    S_2( cutoff+1 : end-cutoff) = 0;

    s_1 = real( ifft( S_1 ) );
    s_2 = real( ifft( S_2 ) );
    
    % avoid log of 0
    s_1 = 20*log10( abs( s_1 ) + 5*eps );
    s_2 = 20*log10( abs( s_2 ) + 5*eps );

    plot( ( taps_to_show - 1302 ) ./ fs .* 1000, s_1, 'k', 'Color', [ .7 .7 .7 ] )
    hold on;
    plot( ( taps_to_show - 1302 ) ./ fs .* 1000, s_2, 'k', 'LineWidth', 2 )
    hold off;

    legend( 'hp', 'lp', 'Location', 'NorthEast' );
    
end
    
ylim( [ -30 0 ] );
xlim( t_to_show .* 1000 );
axis square;
grid on;
xlabel( 't (ms)' );
graph_defaults;


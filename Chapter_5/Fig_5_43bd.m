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

subfigure = 'b'; % for Fig. 5.43(b)
%subfigure = 'd'; % for Fig. 5.43(d)

x_s =   0;
y_s = -.5;
L   = 56;
r_0 = 1.5;
x   = 1;
y   = 0;

no_of_taps = 1024;

fs     = 44100;
c      = 343;
cutoff = 1700; % for Fig. 5.43(d)

% time instances to plot
t_to_show    = [ -0.004 0.01 ]; % in sec
taps_to_show = round( ( t_to_show(1) * fs : t_to_show(2) * fs ) + 1060 );
    
d_alpha_0 = 2 * pi / L;

maximum_delay = ceil( r_0 / c * fs ); % in samples

d = zeros( maximum_delay + 5 * no_of_taps + 2, L );

predelay = round( r_0 / c * fs );

prefilter = wavread( 'wfs_prefilter_100_1800_44100.wav' );

% loop over secondary sources
for l = 1 : L 
    
    alpha_0 = l * d_alpha_0;

    x_0 = r_0 * cos( alpha_0 );
    y_0 = r_0 * sin( alpha_0 );
        
    % angle between virtual source and secondary source 
    alpha_s_prime = atan2( y_0 - y_s, x_0 - x_s );

    % equiv. to Eq. (3.89); warning: this works only if the source is on
    % the y-axis
    if ( alpha_s_prime > 0 )
         continue;
    end
   
    d_xs_x0 = sqrt( ( x_s - x_0 ).^2 + ( y_s - y_0 ).^2 );
 
    % from Eq. (5.14)
    delay = d_xs_x0 / c;

    delay = round( delay * fs );
    
    % time reversal, Sec. 5.6.1
    delay = predelay - delay;
    
    % from Eq. (5.14)
    amplitude = sqrt( r_0 / ( d_xs_x0 + r_0 ) ) / d_xs_x0;
    
    d( delay + 1 : delay + length( prefilter ), l ) = amplitude .* prefilter;
       
end

% put zeros around to have some headroom
d = [ zeros( no_of_taps, L ); d; zeros( no_of_taps, L ) ];

% normalize
d = d ./ max(abs( d( : ) ) );

% initialize s
s = zeros( length( taps_to_show ), 1 );

% this loop evaluates Eq. (5.69)
for l = 1 : L
    
    x_0 = r_0 * cos( (l-1) * d_alpha_0 );
    y_0 = r_0 * sin( (l-1) * d_alpha_0 );
    
    r = sqrt( ( x - x_0 ).^2 + ( y - y_0 ).^2 );
    
    % from Eq. (5.69)
    t = ( taps_to_show./fs - r./c ) .* fs + 1; % in samples

    % Interpolate the impulse responses to find the values at instances
    % t that we are interested in. 
    d_interp = interp1( ( 1 : size( d, 1 ) ), d( :, l ), t, 'linear').';
    
    % from Eq. (5.69)
    s = s + d_interp ./ r;
    
end

figure;

if ( subfigure == 'b' )
    % avoid log of 0
    s = 20*log10( abs( s ) + 5*eps );
    
    plot( ( taps_to_show - 1290 ) ./ fs .* 1000, s, 'k' )
    
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

    plot( ( taps_to_show - 1290 ) ./ fs .* 1000, s_1, 'k', 'Color', [ .7 .7 .7 ] )
    hold on;
    plot( ( taps_to_show - 1290 ) ./ fs .* 1000, s_2, 'k', 'LineWidth', 2 )
    hold off;
    
    legend( 'hp', 'lp', 'Location', 'NorthEast' )
    
end
    
ylim( [ -30 0 ] );
xlim( t_to_show .* 1000 );
axis square;
grid on;
xlabel('t (ms)');
graph_defaults;

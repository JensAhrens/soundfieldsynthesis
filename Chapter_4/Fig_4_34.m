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

subfigure = 'a'; % for Fig. 4.34(a)
%subfigure = 'b'; % for Fig. 4.34(b)

x        = 0;
y        = 1;

y_ref    =  1;
delta_x  = .2;
theta_pw = pi/4;
phi_pw   = pi/2;
fs       = 44100;
c        = 343;
cutoff   = 1800; % for Fig. 4.34(b)

% time instances to plot
t_to_show    = [ -0.004 0.025 ]; % in sec
taps_to_show = round( ( t_to_show(1) * fs : t_to_show(2) * fs ) + 2160 );
    
% create frequency axis
no_of_freq_bins = 2049;
F               = linspace( 0, fs, no_of_freq_bins+1 );
F               = F( 1 : end/2 ).';
F(1)            = F(2); % to avoid numerical instabilities
K               = 2*pi*F / c;

X_0             = -10 : delta_x : 10;
[ x_0, k ]      = meshgrid( X_0, K );
y_0             = 0;
   
% Eq. (A.3a), (A.3b)
k_pw_x = k .* cos( theta_pw ) .* sin( phi_pw );
k_pw_y = k .* sin( theta_pw ) .* sin( phi_pw );

% Eq. (3.79)
D = ( 4*1i .* exp( -1i .* k_pw_y .* y_ref ) ) ./ ...
            besselh( 0, 2, k_pw_y .* y_ref ) .* ...
                    exp( -1i .* k_pw_x .* x_0) ;

d = real( ifft( [ D; conj( flipud( D( 2 : end-1, : ) ) ) ] ) );

% move t = 0 to the center of the buffer
d = circshift( d, [ ( no_of_freq_bins - 1 ) / 2 0 ] );
        
% put zeros around to have some headroom
d = [ zeros( 1024, length( X_0 ) ); d; zeros( 1024, length( X_0 ) ) ];
  
% normalize
d = d ./ max( abs( d( : ) ) );

% initialize s
s = zeros( length( taps_to_show ), 1 );

% this loop evaluates Eq. (5.69)
for index = 1 : length( X_0 )
    
    x_0 = X_0( index );
    
    r = sqrt( ( x - x_0 ).^2 + ( y - y_0 ).^2 );
    
    % from Eq. (5.69)
    t = ( taps_to_show./fs - r./c ) .* fs + 1; % in samples

    % Interpolate the impulse responses to find the values at instances
    % t that we are interested in. 
    d_interp = interp1( ( 1 : size( d, 1 ) ), d( :, index ), t, 'linear').';
    
    % from Eq. (5.67)
    s = s + d_interp ./ r;
    
end

figure;

if ( subfigure == 'a' )
    % avoid log of 0
    s = 20*log10( abs( s ) + 5*eps );
    
    plot( ( taps_to_show - 2150 ) ./ fs .* 1000, s, 'k' )
    
else
    % apply filtering in FFT domain (ugly but serves its purpose...)
    S   = fft( s );
    S_1 = S;
    S_2 = S;
    
    % cutoff frequency in samples
    cutoff = round( cutoff./fs .* size( s, 1 ) ) + 1; % samples

    % set all portions that we want to suppress to zero
    S_1( 1 : cutoff )            = 0; 
    S_1( end-cutoff+1 : end )    = 0;
    S_2( cutoff+1 : end-cutoff ) = 0;

    s_1 = real( ifft( S_1 ) );
    s_2 = real( ifft( S_2 ) );
    
    % avoid log of 0
    s_1 = 20*log10( abs( s_1 ) + 5*eps );
    s_2 = 20*log10( abs( s_2 ) + 5*eps );

    plot( ( taps_to_show - 2150 ) ./ fs .* 1000, s_1, 'k', 'Color', [ .7 .7 .7 ] )
    hold on;
    plot( ( taps_to_show - 2150 ) ./ fs .* 1000, s_2, 'k', 'LineWidth', 2 )
    hold off;

    legend( 'hp', 'lp', 'Location', 'NorthEast' );
    
end
    
ylim( [ -30 0 ] );
xlim( t_to_show .* 1000 );
axis square;
grid on;
xlabel( 't (ms)' );
graph_defaults;

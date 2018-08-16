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

t = 0.013;
v = 150; % m/s
c = 343;
M = v/c;
f = 500; % Hz
omega = 2*pi*f;

alpha_rot = pi/6;

% create spatial grid
X        = linspace( -3, 3, 300 );
Y        = linspace(  0, 6, 300 );
[ x, y ] = meshgrid( X, Y );
z        = 0;

d_ref = 1;
y_0   = 1;
z_0   = 0;

% initialize s
s = zeros( size( x ) );

% this loop evaluates Eq. (5.69) with (5.80)
for x_0 = -3 : 0.1 : 3
    
    disp( [ 'Processing secondary source at x_0 = ' num2str( x_0, '%0.1f' ) ' m.' ] ); 

    r_0 = sqrt( ( x - x_0 ).^2 + ( y - y_0 ).^2 + ( z - z_0 ).^2 );
    tt  = t - r_0./c;
    
    x_s_t = v.*tt;
    
    % Eq. (5.58)
    Delta = sqrt( ( x_0 - x_s_t ).^2 + ( y_0.^2 + z_0.^2 ) .* ( 1 - M.^2 ) );

    % Eq. (5.57)
    tau   = ( M .* ( x_0 - x_s_t ) + Delta ) ./ ( c * ( 1 - M^2 ) );
    
    % from Eq. (5.69) 
    t_tilde     = tt - tau;
    
    x_s_t_tilde = v .* t_tilde;   
    r_tilde     = sqrt( ( x_0 - x_s_t_tilde ).^2 + y_0.^2 + z_0.^2 );
    alpha_tilde = atan2( y_0, x_0 - x_s_t_tilde );
    
    % Eq. (5.60)
    source_field    = 1./(4*pi) .* exp( 1i .* omega .* t_tilde ) ./ Delta;

    % Eq. (5.66)
    source_field_dy = - y_0 ./ Delta .* ( ( 1 - M^2 ) ./ Delta  + (1i*omega)/c ) .* ...
                                                    exp( 1i .* omega .* t_tilde ) ./ Delta;
     
    % Eq. (5.81)
    s_bar_prime     = cos( alpha_tilde - alpha_rot );

    s_bar_prime_dy  =  - ( sin( alpha_tilde - alpha_rot ) .* ( x_0 - x_s_t_tilde ) ) ./ ...
                            ( ( x_0 - x_s_t_tilde ).^2 + y_0.^2 ) .* ...
                                ( 1 - ( ( M .* y_0.^2 ) ./ ( Delta .* ( x_0 - x_s_t_tilde ) ) ) );

    %%%% convolution with sign (Eq. (5.79)) is ommitted for simplicity %%%% 
                                
    % Eq. (5.80)
    d = sqrt( 2*pi*d_ref / ( 1i * omega/c ) ) .* ...
            c / (4*pi) .* ( s_bar_prime  .* source_field_dy + ...
                                          source_field .* s_bar_prime_dy  );
    
    % Eq. (5.69)                              
    s = s + 1/(4*pi) .* d ./ r_0 ;
    
end

% normalize
s = s ./ abs( s( end/2, end/2 ) );

figure;
imagesc( X, Y, real( s ), [ -6 6 ] );
turn_imagesc;
axis square;
colormap gray;

hold on;
plot( -3 : 0.1 : 3, y_0, 'kx' );
hold off;

xlabel( 'x (m)' )
ylabel( 'y (m)' )
graph_defaults;

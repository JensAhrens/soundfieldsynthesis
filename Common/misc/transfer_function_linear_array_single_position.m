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

function [ S ] = transfer_function_linear_array_single_position( x, f )
% this function calculates the transfer function of a dicrete linear
% secondary source distribution to the point x 

y_ref   =  1;
delta_x = .2; 
Nu      = -20 : 20; % spectral repetitions to consider

c = 343;
k = 2*pi*f/c;

theta_pw = pi/4;
phi_pw   = pi/2;

k_pw_x = k * cos( theta_pw ) * sin( phi_pw );
k_pw_y = k * sin( theta_pw ) * sin( phi_pw );

% initialize S
S = zeros( length( f ), 1 );

% this loop evaluates Eq.(4.49)
for nu = Nu

    % argument of D_tilde in Eq.(4.47)
    k_x = (2*pi) / delta_x * nu + k_pw_x;
    
    % initialize G_kx
    G_kx = zeros( length(f), 1 );
    
    % Eq.(C.10)
    G_kx( abs( k_x ) < k ) = -1i/4 * besselh( 0, 2, sqrt( ( k( abs( k_x ) < k ) ).^2 - k_x( abs( k_x ) < k ).^2 ) .* x(2) );
    G_kx( abs( k_x ) > k ) = 1/(2*pi) * besselk( 0, sqrt( k_x( abs( k_x ) > k ).^2 - ( k( abs( k_x ) > k ) ).^2 ) .* x(2) );
    
    % Eq.(4.49)
    S = S + ( 4 * 1i * exp( -1i .* k_pw_y .* y_ref ) ) ./ ...
                ( besselh( 0, 2, k_pw_y .* y_ref ) ) .* ...
                    exp( -1i .* k_x .* x(1) ) .* G_kx;
    
end

end

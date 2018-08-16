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

function [ out ] = Lambda_projection( m, alpha_n, x_0, y_0, k, Q )
% Implements Eq. (5.91)

r_0     = sqrt( x_0^2 + y_0^2 );
alpha_0 = atan2( y_0, x_0 );
out     = zeros( size( k ) );

for l = -Q : Q
    
    out = out + 1i^(-l) .* besselj( l, k.*r_0 ) .* exp( -1i * l * alpha_0 ) .* ...
            ( exp( -1i * alpha_n ) * mathring_w( -1 - m - l, alpha_n ) + ...
              exp(  1i * alpha_n ) * mathring_w(  1 - m - l, alpha_n ) );
    
end

out = -1i/4 .* k .* out;

end

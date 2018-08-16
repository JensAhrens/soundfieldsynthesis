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

function [ out ] = Psi_projection( n, m )
% implements Eq. (E.44)

m = abs( m );

out = pi * sqrt( ( 2*n + 1 ) / (4*pi) * factorial( n - m ) / factorial( n + m ) ) * ...
           ( 2^( -2*m - 1 ) * gamma( 1 + m + n ) ) / ...
          ( gamma( 1/2 + m/2 ) * gamma( 3/2 + m/2 ) * gamma( 1 - m + n ) ) * ...
                 genHyper( [ ( m + n + 1 )/2, ( m - n )/2, m/2 + 1 ], [ m + 1, ( m + 3 )/2 ], 1 );
          
end

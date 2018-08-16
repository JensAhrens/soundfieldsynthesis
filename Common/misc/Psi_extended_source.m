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

function [ out ] = Psi_extended_source( n, m )
% evaluates Eq.(E.36)
%
% Note that this script uses the genHyper function by Ben Barrowes. See
% file Common/genHyper_1.2/README for credits and license.

m = abs( m );

out = pi * ( 2^( -2*m - 1 ) * gamma( 1 + m + n ) ) / ...
          ( gamma( 1/2 + m/2 ) * gamma( 3/2 + m/2 ) * gamma( 1 - m + n ) ) * ...
                 genHyper( [ ( m + n + 1 ) / 2, ( m - n ) / 2, m / 2 + 1 ], [ m + 1, ( m + 3 ) / 2 ], 1 );
          
end

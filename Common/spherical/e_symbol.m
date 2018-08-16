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

function [ out ] = e_symbol( m120, n120 )
%[ out ] = e_symbol( m120, n120 )
%
% E-symbol as defined in Gumerov, Duraiswami, Elsevier 2004, eq.(3.2.27)
% 
%       / m_1 m_2 m \
%   E  |             |
%       \ n_1 n_2 n /
%
% m120: [ m_1 m_2 m ]
% n120: [ n_1 n_2 n ]
%
% Author: Jens Ahrens

out = 4 * pi * epsilon( m120(1) ) * epsilon( m120(2) ) * epsilon( m120(3) ) * ...
            Wigner3j( n120, [ 0 0 0 ] ) * Wigner3j( n120, m120 );

end

function [ out ] = epsilon( m ) 

out = 1;

if ( m > 0 )
    out = (-1)^m;
end

end
function [ Ynm ] = sphharm( n, m, beta, alpha, type )
% spherical harmonics of degree n and order m.
%
% [ Ynm ] = sphharm( n, m, beta, alpha, type );
%
% n           - spherical harmonic degree
% m           - spherical harmonic order
% beta        - colatitude to be calculated
% alpha       - azimuth to be calculated 
% type        - 'complex' (default) or 'real'
%
% alpha and beta can be arrays but have to be of same size or one of them
% has to be a scalar.
%
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

if ( nargin < 5 )
    type = 'complex';
end

if ( n < 0 )
    error( 'Degree(n) must not be negative.' )
end

if ( n < abs(m) )
    warning( 'Absolute value of order m must be less than or equal to the degree n.' ); 
    Ynm = zeros( size( alpha ) );
    return;
end

Lnm = asslegendre( n, abs(m), cos( beta ) );

factor_1 = ( 2*n + 1 ) / ( 4*pi );
factor_2 = factorial( n - abs(m) ) ./ factorial( n + abs(m) );

% complex spherical harmonics
if ( strcmp( type, 'complex' ) )
    
    Ynm = (-1).^m .* sqrt( factor_1 .* factor_2 ) .* Lnm .* exp( 1i .* m .* alpha );

% real valued spherical harmonics
elseif( strcmp( type, 'real' ) )
        
    if ( m ~= 0 )
        factor_1 = 2 * factor_1;
    end
    
    if ( m < 0 )
        Ynm = (-1).^m .* sqrt( factor_1 .* factor_2 ) .* Lnm .* sin( m .* alpha );
    else % m >= 0
        Ynm = (-1).^m .* sqrt( factor_1 .* factor_2 ) .* Lnm .* cos( m .* alpha );
    end
    
else
    error( 'Unknown type.' );    
    
end

end

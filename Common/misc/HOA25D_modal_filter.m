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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The author of this script is Sascha Spors. It is also part of the Sound
% Field Synthesis Toolbox, which can be obtained from
% http://dev.qu.tu-berlin.de/projects/sfs.

function [ hd ] = HOA25D_modal_filter( src_type, R, order, fs, c, r_s )
% src_type: either 'pw' or 'sw' for plane or spherical waves
% r_s is the source distance for if src_type = 'sw'; not needed for 'pw';

% compute normalized roots/zeros of spherical Hankel function
B = zeros( 1, order + 2 );
A = B;

for n = 0 : order
    B( n+1 ) = factorial( 2*order - n ) / ( factorial( order - n ) * factorial( n ) * 2^( order - n ) );
end

B = fliplr( B );
% find zeros/roots
z = roots( B );

% compute SOS coefficients of modal driving function
if ( strcmp( src_type, 'pw' ) )
    A(2) = 1;
    p = roots( A );
    [ sos, g ] = zp2sos( p, z*c/R, 2, 'down', 'none' );
elseif ( strcmp( src_type, 'sw' ) )
    A(1) = 1;
    [ sos, g ] = zp2sos( z * c/r_s, z * c/R, 1, 'up', 'none' ); 
end

% transform coefficients
for n = 1 : size( sos, 1 )
    
    [ bz, az ] = bilinear( sos( n, 1 : 3 ), sos( n, 4 : 6 ), fs, 1000 );
    
    if ( length( bz ) == 2 )
        sos( n, 2 : 3 ) = bz;
        sos( n, 5 : 6 ) = az;
    else
        sos( n, 1 : 3 ) = bz;
        sos( n, 4 : 6 ) = az;
    end
end

% realize FOS/SOS as DF-II structure
hd = dfilt.df2sos( sos ); 

if ( strcmp( src_type, 'pw' ) )
    hd.ScaleValues( end ) = 2*(-1)^order;

elseif ( strcmp( src_type, 'sw' ) )
    hd.ScaleValues( end ) = 1/(2*pi*R);

end

end

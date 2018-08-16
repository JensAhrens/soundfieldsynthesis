function [ out, weight ] = normalize_audio_signal( in )

v = max( abs( in( : ) ) );

weight = 0.99 / v( 1 );

out = in * weight;

end
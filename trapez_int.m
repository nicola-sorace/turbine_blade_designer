% Integration using trapezoidal method
% TODO better method?

function A = trapez_int(ps, vs)
    A = sum((ps(2:end)-ps(1:end-1)) .* (vs(1:end-1)+vs(2:end))/2 );
end
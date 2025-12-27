function newAlphaNR = emdlab_flib_updateAlphaNR(alphaNR, relativeError)
if length(relativeError)>2
    if relativeError(end) > 0.8*relativeError(end-1)
        newAlphaNR = max(0.95*alphaNR-2e-3*rand,0.5);
    else
        newAlphaNR = min(1.05*alphaNR+2e-3*rand,0.9);
    end
else
    newAlphaNR = alphaNR;
end
end
function out = randOoM(range,stream,baseExp)
    out = 10.^(range*(rand(stream) - .5 + baseExp));
end
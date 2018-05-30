function ttl = stimttlfun(ontime)
    if ontime > 10
        ttl = 0;
    else
        ttl = ceil(ontime/1.35);
end

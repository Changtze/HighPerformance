function minDist = min_particle_distance(data)
    uniqueP = unique(data.pNum);
    minDist = inf;
    
    for i = 1:length(uniqueP)
        for j = i+1:length(uniqueP)
            p1 = data(data.pNum == uniqueP(i), :);
            p2 = data(data.pNum == uniqueP(j), :);
            
            for k = 1:height(p1)
                for m = 1:height(p2)
                    dist = sqrt((p1.x(k) - p2.x(m))^2 + (p1.y(k) - p2.y(m))^2 + (p1.z(k) - p2.z(m))^2);
                    minDist = min(minDist, dist);
                end
            end
        end
    end
end

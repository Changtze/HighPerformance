function plot_particle_distance(data, p1, p2)
    % Ensure variable names are valid identifiers
    data.Properties.VariableNames = matlab.lang.makeValidName(data.Properties.VariableNames);
    
    % Find common time points
    t_common = intersect(data.t(data.pNum == p1), data.t(data.pNum == p2));
    distances = zeros(size(t_common));

    for i = 1:length(t_common)
        t = t_common(i);
        pos1 = data{data.pNum == p1 & data.t == t, {'x', 'y', 'z'}};
        pos2 = data{data.pNum == p2 & data.t == t, {'x', 'y', 'z'}};

        if ~isempty(pos1) && ~isempty(pos2)
            distances(i) = sqrt(sum((pos1 - pos2).^2));
        end
    end

    plot(t_common, distances, '-o');
    xlabel('Time');
    ylabel('Distance');
    title(sprintf('Distance Between Particle %d and %d Over Time', p1, p2));
end

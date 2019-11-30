function hitting_set = greedy_hitting_set(sets)
set_elements = unique(cell2mat(cellfun(@(x)x(:).',sets, 'UniformOutput', false)));

remaining_sets = sets;
remaining_elements = set_elements;
hitting_set = [];
while ~isempty(remaining_sets)
    greedy_element = [];
    greedy_intersections = [];
    for i = 1:length(remaining_elements)
        intersections = cellfun(@(x) intersect(x,remaining_elements(i)),...
            remaining_sets,'UniformOutput',false);
        intersections = cell2mat(cellfun(@(x) ~isempty(x),intersections,'UniformOutput',false));
        if sum(intersections) > sum(greedy_intersections)
            greedy_element = remaining_elements(i);
            greedy_intersections = intersections;
        end
        
    end
    hitting_set = [hitting_set, greedy_element];
    remaining_elements = remaining_elements(remaining_elements ~= greedy_element);
    remaining_sets = remaining_sets(~greedy_intersections);
end

hitting_set = sort(hitting_set);

end


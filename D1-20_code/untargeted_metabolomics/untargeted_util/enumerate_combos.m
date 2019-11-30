function [donor_combos,compound_combos] = enumerate_combos(donors,compounds)

donor_combos = [];
compound_combos = {};
for i = 1:length(donors)
    for j = 1:length(compounds)
        if ~strcmp(compounds{j},'DMSO')
            donor_combos = [donor_combos;donors(i)];
            compound_combos = [compound_combos;compounds{j}];
        end
    end
end

end


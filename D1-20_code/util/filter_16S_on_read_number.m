function filtered_manifest = filter_16S_on_read_number(manifest,read_cutoff)
%This script filters out samples with less than a certain number of reads
num_reads = nan(size(manifest,1),1);
for i = 1:size(manifest,1)
    num_reads(i) = sum(manifest.asv{i}.Variables);
end
filtered_manifest = manifest(~(num_reads < read_cutoff), :);

end


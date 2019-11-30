function manifest = wrapped_16S_import()

asv_file = '16S_data/MDM_2019_16S_HTP_asv_table.txt';
order_file = '16S_data/MDM_2019_16S_HTP_table_4.txt';
family_file = '16S_data/MDM_2019_16S_HTP_table_5.txt';
genus_file = '16S_data/MDM_2019_16S_HTP_table_6.txt';
species_file = '16S_data/MDM_2019_16S_HTP_table_7.txt';
manifest_file = '16S_data/HTP_16S_analysis_manifest.xlsx';
primary_biomass_file = '16S_data/primary_biomass_experiment.csv';
secondary_biomass_file = '16S_data/secondary_biomass_experiment.csv';

manifest = import_16S_results(asv_file,species_file,genus_file,...
    family_file,order_file,manifest_file,primary_biomass_file,secondary_biomass_file);


end


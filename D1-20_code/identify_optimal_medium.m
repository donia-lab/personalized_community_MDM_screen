%This script identifies the best media according to a number of different
%metrics

%% Import data and compute noise models
clear;clc

manifest = wrapped_16S_import();

filtered_manifest = filter_16S_on_read_number(manifest,1e4);

%Compute noise models for ENDS

donor_manifest_loc = 'QTOF_data/20x20_donor_manifest.xlsx';
donor_loc = 'QTOF_data/donor_data.csv';
istd = 'voriconazole';
plt = 0;
noise_models = estimate_noise_models(donor_manifest_loc,donor_loc,istd,plt);

%ENDS parameters and storage vector
alpha = 0.01;
n = 3;
r = 1;
only_ENDS = false;
read_cutoff = 0;

%% Compute 16S scores
score_manifest = ...
    compute_16S_scores(filtered_manifest,r,alpha,n,noise_models,read_cutoff,only_ENDS);

%% Get highest scoring media for different metrics

media_manifest = score_manifest(~contains(score_manifest.media,{'feces','none'}),:);

[optimal_ENDS_media, optimal_ENDS_value] = find_highest_medium(media_manifest,'ENDS');
[optimal_richness_media, optimal_richness_value] = find_highest_medium(media_manifest,'richness');
[optimal_shannon_media, optimal_shannon_value] = find_highest_medium(media_manifest,'shannon');
[optimal_JS_media, optimal_JS_value] = find_highest_medium(media_manifest,'JS');
[optimal_shared_1_pct_media, optimal_shared_1_pct_value] = find_highest_medium(media_manifest,'shared_asv_above_1_pct');



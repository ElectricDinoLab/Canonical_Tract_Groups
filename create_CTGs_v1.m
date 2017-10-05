
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code creates Canonical Tract Groups (CTGs) from functional and
% structural data, and then summarizes the connectivity (F & S) values
% within those CTGs.
%
% Current version uses the JHU Tract atlas, and a sub-parcellated version
% of the Harvard-Oxford atlas with 471 cortical ROIs.
%
% 1st version:  Simon W Davis, Amanda Syzmanski
%               Duke University
%
%%%%%%%%%%
% First step: assess the overlap between tracts and ROIs within a given atlas.
%%%%%%%%%%

clear all
compute_binaries = 1;
compute_TGs = 1;

if compute_binaries == 1
    
    tract_atlas_file = sprintf('/Volumes/Cabeza/MemEX.01/Scripts/JHU-ICBM-tracts-maxprob-thr0-2mm.nii');
    atlas_file = sprintf('/Volumes/Cabeza/MemEX.01/Scripts/HOAsp.nii');
    
    % change to 100 or 471
    tract_atlas471 = spm_read_vols(spm_vol(tract_atlas_file));
    atlas = spm_read_vols(spm_vol(atlas_file));
    
    overlap = zeros(20,471); % change to 100 or 471 based on atlas used
    
    for i = 1:20 % run through all tract groups (left and right)
        for j = 1:471 %change to 100 or 471
            tract_atlas_mod = tract_atlas471;
            tract_atlas_mod(tract_atlas_mod ~= i) = 0;
            
            atlas_mod = atlas;
            atlas_mod(atlas_mod ~= j) = 0;
            
            count = nnz(atlas_mod .* tract_atlas_mod);
            overlap(i, j) = count;
            
        end
    end
    
    % Specify the amount of overlap between the tract and the ROI. Greater
    % value for n will give a more conservative map.
    n = 10;
    count_binary = overlap;
    count_binary(count_binary < n) = 0;
    count_binary(count_binary >= n) = 1;
end

if compute_TGs == 1
    tract_atlas471 = zeros(471,471, 20); % change to 100 or 471 based on atlas used for a given tract
    for k = 1:size(count_binary, 1)
        
        % which regions have a nonzero value / what column/y values have a nonzero
        % value / determine the y indices that have nonzero value
        x = count_binary(k, :);
        test = find(x>0);
        % generate all possible pairs of regions
        for blah = 1:size(test, 2)
            for blah2 = 1:size(test, 2)
                tract_atlas471(test(blah), test(blah2), k) = 1;
            end
            save /Volumes/Cabeza/MemEX.01/Scripts/tract_groups_mar17_471_atlas.mat tract_atlas471   % change to 100 or 471
        end
        
    end    
end


%%%%%%%%%%%%%%%%%%
% Second step: Compute means for a given CTG, for a given modality
%%%%%%%%%%%%%%%%%%

subjects = {'35014','35106','35136','35148','35172','35186','35217','35228','35271','35275','35310','35311','35312','35320','35336','35350','35389','35396','35425','35446','35450','35453','35486','35493','35505','35508','35552','35576','35581','35583','35609','35611','35653','35655','35659','35714','35793','35800','35808','35819','35823','35832','35923','35976','35634','35649','35663','35686','35687','35707','35721','35722','35729','35849','35856','35860','35900','35917','35950','35647','35669','35967'};
ages = {'O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','O','Y','Y','Y','Y','Y','Y','Y','Y','Y','Y','Y','Y','Y','Y','Y','Y','Y','Y'};

% flags to turn on/off functions.
TGs_DTI = 0; % tract groups, structure
TGs_Func = 1; % tract groups, hits

% The DTI CTGs
if TGs_DTI == 1  
    load tract_groups_mar17_471_atlas.mat
    for i = 1:20        % "i" is the counter for the tract group (UF, ILF, etc)
        parfor j = 1:size(subjects, 2)
            try
                sub_func_file = sprintf('/Volumes/Cabeza/MemEX.01/Data/Preprocessing/%s/%s_FA471.csv', subjects{j}, subjects{j}); % if you want SIFTed tracts, add "_SIFT" to the end of the csv
                subject_connectome = dlmread(sub_func_file);                    % load subject data
                tract_groups_filter = squeeze(tract_atlas471(:,:,i));           
                filtered = triu(tract_groups_filter) .* subject_connectome;     % extract the upper triangle of the CTG filter and multiply by the subjects' data
                a = reshape(filtered,(size(tract_groups_filter,1)*size(tract_groups_filter,1)),1);
                TG_DTI_471(j,i) = nanmean(a(a~=0));                             % save the sum of all cells to file
            catch  
                fprintf('DTI error on TG %d subject %s\n', i, subjects{j});
            end
        end
        fprintf('Finished with TG %d\n', i);
    end
    fprintf('Finished with DTI\n');
end

% The Functional CTGs (one block item memory, one block source memory)
if TGs_Func == 1  
    load tract_groups_mar17_471_atlas.mat
    for i = 1:20        % "i" is the counter for the tract group (UF, ILF, etc)
        for j = 1:size(subjects, 2)
            try
                sub_func_file = sprintf('/Volumes/Cabeza/MemEX.01/Analysis/FunctionalConnectomes/HOA471/%s_item_RSA_471matrix.csv', subjects{j});
                subject_connectome(j,:,:) = dlmread(sub_func_file);             % load subject data
                tract_groups_filter = squeeze(tract_atlas471(:,:,i));
                filtered = triu(tract_groups_filter) .* subject_connectome;     % extract the upper triangle of the CTG filter and multiply by the subjects' data
                b = reshape(filtered,(size(tract_groups_filter,1)*size(tract_groups_filter,1)),1);
                TG_fMRI_471_item(j,i) = nanmean(b(b~=0));                       % save the sum of all cells to file
            catch
                fprintf('Item error on TG %d subject %s\n', i, subjects{j});
            end
        end
        fprintf('Finished with TG %d\n', i);
    end
    fprintf('Finished with Item');
    clear j i sub_func_file subject_connectome tract_groups_filter
    
    load tract_groups_mar17_471_atlas.mat
    for i = 1:20        % "i" is the counter for the tract group (UF, ILF, etc)
        parfor j = 1:size(subjects, 2)
            try
                sub_func_file = sprintf('/Volumes/Cabeza/MemEX.01/Analysis/FunctionalConnectomes/HOA471/%s_source_RSA_471matrix.csv', subjects{j});
                subject_connectome = dlmread(sub_func_file);                    % load subject data
                tract_groups_filter = squeeze(tract_atlas471(:,:,i));
                filtered = triu(tract_groups_filter) .* subject_connectome;     % extract the upper triangle of the CTG filter and multiply by the subjects' data
                c = reshape(filtered,(size(tract_groups_filter,1)*size(tract_groups_filter,1)),1);
                TG_fMRI_471_source(j,i) = nanmean(c(c~=0));                     % save the sum of all cells to file
            catch
                fprintf('Source error on TG %d subject %s\n', i, subjects{j});
            end
        end
        fprintf('Finished with TG %d\n', i);
    end
    fprintf('Finished with Source\n');
    clear j i sub_func_file subject_connectome tract_groups_filter
    
    % Add in some relevant behavioral details to keep track!
    subject.IDs = subjects';
    subject.AGES = ages';
    
    save /Volumes/Cabeza/MemEX.01/Analysis/TGs_471_atlas_Apr16_FORNIX.mat TG_fMRI_471_source TG_fMRI_471_item TG_DTI_471 subject
end

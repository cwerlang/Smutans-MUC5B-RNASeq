%% Differential expression analysis for data from Werlang et. al. Nat. Micro. (2021)
% DOI: 10.1038/s41564-021-00876-1
% Original method/code from Anders, S. & Huber, W.  Genome Biol. (2010)
% Adapted by C. Werlang and K. Wheeler

%%
clear 
clc
close all

%% Load in working directory & data output from R code
% In MATLAB, read in the expression dataset using the readtable function.
% Columns for exprsMat are: gene, sample1_counts, sample2_counts,...

exprsMat = readtable(fullfile(pwd,'data','interim',...
    'raw_read_counts.csv'),'ReadRowNames',true);

% Inspect the first 10 rows in the table with the counts for all samples
exprsMat(1:10, :)

%% Load in file with structure of data

replicates = readtable(fullfile(pwd,...
    'experiment_design_replicate_matrix.csv'),'ReadRowNames',true)
comparisons = readtable(fullfile(pwd,...
    'experiment_design_comparison_matrix.csv'),'ReadRowNames',true)

% set up empty array to write data labels and filenames into later
%output_list = zeros(sum(sum(comparisons{:,:})),2);
output_list = {};
counter = 0;
%% Iterate over the data to generate different comparisons
for i = 1:size(comparisons,1)
    for j = 1:size(comparisons,2)
        close all
        if comparisons{i,j} == 0 % don't compare samples that don't want to be compared
            continue
        end
        % now we are comparing samples
        Mediax = logical(replicates{i,:}); % set the row as baseline
        Treatx = logical(replicates{j,:}); % set the column as treatment
        % do some housekeeping to make sure the script knows which samples are which
        whole = logical(Mediax+Treatx);
        relevant_replicates = [find(Mediax) find(Treatx)]; % find nonzero elements in logical masks
        Media_mask = Mediax(relevant_replicates);
        Treat_mask = Treatx(relevant_replicates);
        % set the file name 
        NAME = strcat(replicates.Properties.RowNames{j},'vs',replicates.Properties.RowNames{i});
        FILE_NAME = strcat('RelExp_',NAME,'.csv');
        FILE_OUTPUT = fullfile(pwd,'data','interim',FILE_NAME);
        % write the file name into matrix for datakeeping
        counter = counter+1;
        output_list{counter,1} = NAME;
        output_list{counter,2} = FILE_NAME;
        
        %% Estimating Library Size Factor
        % The expectation values of all gene counts from a sample are proportional
        % to the sample's library size. The effective library size can be estimated
        % from the count data.

        %%
        % Compute the geometric mean of the gene counts (rows in |exprsMat|) across
        % all samples in the experiment as a pseudo-reference sample.

        pseudo_ref_sample = geomean(exprsMat{:,whole},2);

        %%
        % Each library size parameter is computed as the median of the ratio of
        % the sample's counts to those of the pseudo-reference sample.

        nzi = pseudo_ref_sample>100; % ignore genes with a geometric mean less than 100
        ratios = bsxfun(@rdivide, exprsMat{nzi,whole}, pseudo_ref_sample(nzi));
        sizeFactors = median(ratios, 1);

        %% 
        % The counts can be transformed to a common scale using size factor
        % adjustment.

        base_exprsMat = exprsMat;
        base_exprsMat{:,whole} = bsxfun(@rdivide,exprsMat{:,whole},sizeFactors);

        %%
        % Use the |boxplot| function to inspect the count distribution of the
        % pgm-treated and lb-treated samples and the size factor adjustment.

        figure
        subplot(2,1,1)
        maboxplot(log2(exprsMat{:,whole}), 'title','Raw Read Counts',...
                                          'orientation', 'horizontal') 
        subplot(2,1,2)
        maboxplot(log2(base_exprsMat{:,whole}), 'title','Size Factor Adjusted Read Counts',...
                                               'orientation', 'horizontal')

        %% Estimate the gene abundance 
        % To estimate the gene abundance for each experimental condition you use 
        % the average of the counts from the samples transformed to the common 
        % scale. (Eq. 6 in [6])

        mean_Media = mean(base_exprsMat{:,Mediax}, 2); 
        mean_Treat = mean(base_exprsMat{:,Treatx}, 2); 

        % remove infinite and Nan values -- replace with 0
        mean_Treat(isinf(mean_Treat)|isnan(mean_Treat)) = 0;
        mean_Media(isinf(mean_Media)|isnan(mean_Media)) = 0;

        % consider the dispersion
        disp_Media = std(base_exprsMat{:,Mediax},0,2) ./ mean_Media;
        disp_Treat = std(base_exprsMat{:,Treatx},0,2) ./ mean_Treat;

        %% plot mean vs dispersion on a log-log scale

        figure
        loglog(mean_Treat,disp_Treat,'r.');
        hold on;
        loglog(mean_Media,disp_Media,'b.');
        xlabel('log2(Mean)');
        ylabel('log2(Dispersion)');
        legend(NAME,'Media','Location','southwest');


        %%
        % Plot the log2 fold changes against the base means using the |mairplot|
        % function. A quick exploration reflects ~10 differentially expressed genes
        % (20 fold change or more), though not all of these are significant due to
        % the low number of counts compared to the sample variance.

        mairplot(mean_Treat(nzi),mean_Media(nzi),'Labels',exprsMat.Properties.RowNames,'Factor',5)

        %%
        % MA Plot
        mairplot(mean_Treat(nzi),mean_Media(nzi),'Labels',exprsMat.Properties.RowNames,'Type','MA');
        set(get(gca,'Xlabel'),'String','mean of normalized counts')
        set(get(gca,'Ylabel'),'String','log2(fold change)')

        %% Estimating Negative Binomial Distribution Parameters
        % In the model, the variances of the counts of a gene are considered as the
        % sum of a shot noise term and a raw variance term. The shot noise term is
        % the mean counts of the gene, while the raw variance can be predicted from
        % the mean, i.e., genes with a similar expression level have similar
        % variance across the replicates (samples of the same biological
        % condition). A smooth function that models the dependence of the raw
        % variance on the mean is obtained by fitting the sample mean and variance
        % within replicates for each gene using local regression function. 

        %%
        % Compute sample variances transformed to the common scale for pgm-treated
        % samples. (Eq. 7 in [6]) 

        var_Treat = var(base_exprsMat{:,Treatx}, 0, 2);
        var_Media = var(base_exprsMat{:,Mediax}, 0, 2);


        %%
        % Estimate the shot noise term. (Eq. 8 in [6])
        z1 = mean_Treat * mean(1./sizeFactors(Treat_mask));
        z2 = mean_Media * mean(1./sizeFactors(Media_mask));


        % The helper function |estimateNBVarFunc| returns an anonymous function
        % that maps the mean estimate to an unbiased raw variance estimate. Bias
        % adjustment due to shot noise and multiple replicates is considered in the
        % anonymous function.  

        raw_var_func_Media = estimateNBVarFunc(mean_Media,var_Media,sizeFactors(Media_mask));
        raw_var_func_Treat = estimateNBVarFunc(mean_Treat,var_Treat,sizeFactors(Treat_mask));

        %%
        % Use the anonymous function |raw_var_func_A| to calculate the sample
        % variance by adding the shot noise bias term to the raw variance. (Eq.9 in
        % [6])  
        var_fit_Treat = raw_var_func_Treat(mean_Treat) + z1;
        var_fit_Media = raw_var_func_Media(mean_Media) + z2;

        %%
        % Plot the sample variance to its regressed value to check the fit of the
        % variance function. 
        figure
        loglog(mean_Media, var_Media, '*')
        hold on
        loglog(mean_Media, var_fit_Media, '.r')
        ylabel('Base Variances')
        xlabel('Base Means')
        title(strcat('Dependence of the Variance on the Mean for untreated Samples'))


        %%
        % Plot the sample variance to its regressed value to check the fit of the
        % variance function. 
        figure
        loglog(mean_Treat, var_Treat, '*')
        hold on
        loglog(mean_Treat, var_fit_Treat, '.r')
        ylabel('Base Variances')
        xlabel('Base Means')
        title(strcat('Dependence of the Variance on the Mean for', NAME ,'-Treated Samples'))

        %%
        % The fit (red line) follows the single-gene estimates well, even though
        % the spread of the latter is considerable, as one would expect, given that
        % each raw variance value is estimated from only three values (three
        % pgm-treated replicates).

        %% Empirical Cumulative Distribution Functions
        % As RNA-seq experiments typically have few replicates, the single-gene
        % estimate of the base variance can deviate wildly from the fitted value.
        % To see whether this might be too wild, the cumulative probability for the
        % ratio of single-gene estimate of the base variance to the fitted value is
        % calculated from the chi-square distribution, as explained in reference
        % [6].

        % %%
        % % Compute the cumulative probabilities of the variance ratios of
        % % treated samples.
        % degrees_of_freedom = sum(Treatx) - 1;
        % var_ratio = var_Treat ./ var_fit_Treat;
        % pchisq = chi2cdf(degrees_of_freedom * var_ratio, degrees_of_freedom);
        %%
        % Compute the cumulative probabilities of the variance ratios of
        % untreated samples.
        degrees_of_freedom = sum(Mediax) - 1;
        var_ratio = var_Media ./ var_fit_Media;
        pchisq = chi2cdf(degrees_of_freedom * var_ratio, degrees_of_freedom);

        %%
        % Compute the empirical cumulative density functions (ECDF) stratified by
        % base count levels, and show the ECDFs curves. Group the counts into seven
        % levels.
        count_levels = [0 3 12 30 65 130 310];
        labels = {'0-3','4-12','13-30','31-65','66-130','131-310','> 311'};
        %grps = sum(bsxfun(@ge,mean_Treat,count_levels),2); % stratification
        grps = sum(bsxfun(@ge,mean_Media,count_levels),2); % stratification

        %%
        % The ECDF curves of count levels greater than 3 and below 130 follows the
        % diagonal well (black line). If the ECDF curves are below the black line,
        % variance is underestimated. If the ECDF curves are above the black line,
        % variance is overestimated [6]. For very low counts (below 3), the
        % deviations become stronger, but at these levels, shot noise dominates.
        % For the high count cases, the variance is overestimated. The reason might
        % be there are not enough genes with high counts. Get the number of genes
        % in each of the count levels.

        array2table(accumarray(grps,1),'VariableNames',{'Counts'},'RowNames',labels)

        %%
        % Increasing the sequence depth, which in turn increases the number of
        % genes with higher counts, improves the variance estimation.

        %% Testing for Differential Expression
        % Having estimated and verified the mean-variance dependence, you can test
        % for differentially expressed genes between the samples from the PGM- and
        % LB- treated conditions. Define, as test statistic, the total counts in
        % each condition, |k_Treat| and |k_Media|:

        k_Media = sum(exprsMat{:, Mediax}, 2);
        k_Treat = sum(exprsMat{:, Treatx}, 2); 

        %%
        % Parameters of the new negative binomial distributions for count sums
        % |k_Treat| can be calculated by Eqs. 12-14 in [6]: 

        pooled_mean = mean(exprsMat{:,whole},2);

        mean_k_Treat = pooled_mean * sum(sizeFactors(Treat_mask));
        var_k_Treat = mean_k_Treat + raw_var_func_Treat(pooled_mean) * sum(sizeFactors(Treat_mask).^2);

        mean_k_Media = pooled_mean * sum(sizeFactors(Media_mask));
        var_k_Media = mean_k_Media + raw_var_func_Media(pooled_mean) * sum(sizeFactors(Media_mask).^2);



        %%
        % Compute the p-values for the statistical significance of the change from
        % PGM-treated condition to LB-treated condition. The helper function
        % |computePVal| implements the numerical computation of the p-values
        % presented in the reference [6].

        res = table(exprsMat.Properties.RowNames,'VariableNames',{'Gene'});
        res.pvals = computePVal(k_Treat, mean_k_Treat, var_k_Treat, k_Media, mean_k_Media, var_k_Media);

        %%
        % You can empirically adjust the p-values from the multiple tests for false
        % discovery rate (FDR) with the Benjamini-Hochberg procedure [7] using the
        % |mafdr| function.

        res.p_fdr = mafdr(res.pvals, 'BHFDR', true);

        %%
        % Determine the fold change estimated from the LB-treated to the
        % PGM-treated condition.

        fold_change = mean_Treat ./ mean_Media;

        %% 
        % Determine the base 2 logarithm of the fold change.

        res.log2_fold_change = log2(fold_change);

        %%
        figure
        xval = 1:length(fold_change);
        %xval is number of genes (rows - 1 on excel sheet)
        scatter(xval,res.log2_fold_change,12,'filled','k')

        moreThanThreshold = res.log2_fold_change > log2(2); % Logical indexes.
        lessThanThreshold = res.log2_fold_change < log2(1/2);
        % Extract those over the threshold into new arrays.
        over_x = xval(moreThanThreshold);
        over_y = res.log2_fold_change(moreThanThreshold);

        under_x = xval(lessThanThreshold); 
        under_y = res.log2_fold_change(lessThanThreshold);
        % Now plot those over 100 with red stars over the first set of points.
        hold on; % Important so you don't blow the first plot away.
        scatter(over_x, over_y,20,'filled','r');
        hold on; % Important so you don't blow the first plot away.
        scatter(under_x, under_y,20,'filled','b');

        % plot lines for up and downregulated
        yline(1,'--r');
        yline(-1,'--b');

        xlabel('Genome Locus', 'FontSize', 20);
        ylabel('Log_{2}(Fold change in RNA)', 'FontSize', 20);
        set(gca, 'FontSize',20);
        title(strcat(NAME))

        %%
        % Volcano Plot
        mavolcanoplot(log2(mean_Treat),log2(mean_Media),res.p_fdr,'Labels',exprsMat.Properties.RowNames)

        %%
        % Create table with statistics about each gene
        %geneTable = table(exprsMat.Properties.RowNames,fold_change,res.p_fdr, 'VariableNames',{'Gene','FC','fdr'});
        %geneTable.Properties.RowNames = exprsMat.Properties.RowNames;
        writetable(res,FILE_OUTPUT)

        %%
        % You can identify up- or down- regulated genes (2 fold change and fdr < .05).
        up_idx = find(res.log2_fold_change >= 1 & res.p_fdr < .1);
        numel(up_idx)
        %%
        down_idx = find(res.log2_fold_change <= -1 & res.p_fdr < .1);
        numel(down_idx)
        
    end
end

%% write list of generated files into a csv
writecell(output_list, fullfile(pwd, 'data','interim','file_list_relative_expression.csv'));
%% References
% [1] Li, H., et al., "Determination of Tag Density Required for Digital
%     Transcriptome Analysis: Application to an Androgen-Sensitive Prostate
%     Cancer Model", PNAS, 105(51):20179-84, 2008.
%
% [2] Langmead, B., Trapnell, C., Pop, M., and Salzberg, S.L., "Ultrafast
%     and Memory-efficient Alignment of Short DNA Sequences to the Human
%     Genome", Genome Biology, 10(3):R25, 2009.
%
% [3] Li, H., et al., "The Sequence Alignment/map (SAM) Format and
%     SAMtools", Bioinformatics, 25(16):2078-9, 2009.
%
% [4] Mortazavi, A., et al., "Mapping and quantifying mammalian
%     transcriptomes by RNA-Seq", Nature Methods, 5:621-8, 2008.
%
% [5] Robinson, M.D. and Oshlack, A., "A Scaling Normalization method for
%     differential Expression Analysis of RNA-seq Data", Genome Biology,
%     11(3):R25, 2010.
%
% [6] Anders, S. and Huber, W., "Differential Expression Analysis for
%     Sequence Count Data", Genome Biology, 11(10):R106, 2010.
%
% [7] Benjamini, Y. and Hochberg, Y., "Controlling the false discovery
%     rate: a  practical and powerful approach to multiple testing",
%     Journal of the Royal Statistical Society, 57(1):289-300, 1995.
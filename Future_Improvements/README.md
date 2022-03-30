# Improvements

This README file describes the improvments that is suggested by the smogn_testing_version.ipynb file. it will be incorporated to the main tool after testing.

## Changes:

under “The core algorithm of active learning”:
- Optional data preprocessing options include SMOTE oversampling for large sets of unevenly distributed data.
- The type of the XGBoost algorithm to optimize is chosen (options: „gbtree“ and „dart“). 

Under “Performance analysis using cross-validation”:

The performance of two algorithm options on original and SMOTE oversampled data was compared using the same k-fold cross-validation setup on the efficiency optimization dataset (1017 points). For SMOGN comparison the entire dataset was used for k-fold crossvalidation analysis, while the comparison between “gbtree” and “dart” algorithms occurred additionally a smaller dataset (k-fold on 360 datapoits, rest – for validation). For the “dart” algorithm additional parameters were included into the parameter grid for optimization: rate_drop [0.02, 0.2], skip_drop [0.01, 0.1]. The resulting mean absolute error, r2 and spearman correlation distributions are presented in Supplementary Figures 1 and 2.

Under Supplementary Note 2: 
- add DART vs default (‘gbtree’) option
- add SMOGN data preprocessing (smogn package [3])


## Explanations:
Active learning in biological applications is further known to suffer from overfitting that results from the small number of labeled samples (cite [2]). To address this issue we included the XGBoost implementation of DART algorithm to our workflow as an option. This algorithm introduces dropout – a widely used regularization method from the field of deep neural networks – boosted trees (cite [1]). 
To additionally correct for the uneven data distribution that can occur at late stages of optimization, optional Synthetic Minority Over-Sampling Technique for Regression with Gaussian Noise (SMOGN, [3]) input preprocessing was included into our workflow. This algorithm allows for synthetic oversampling of underrepresented datapoints based on the original dataset. 
As expected, the use of the DART algorithm did not significantly improve predictions for models trained on the large dataset (data not shown). However, in case of the small dataset optimal parameters allowed for lowering of the prediction variance between cross-validation folds as well as MAE, thus suggesting overfitting reduction (see Supplemental Figure 1). It is worth noting that the DART model training is slower due to the missing use of the prediction buffer and the optimization parameters need to be adapted to the sample sizes. Thus, the overfitting optimization remains a question of gain versus time ratio for each experimental setup.
As can be seen from the presented figures (figure 2) additional data preprocessing using SMOGN significantly improved the performance of the respective models within all tested metrics especially in the high-yield range. This can be explained by the fact that the oversampling of this data resulted in model bias towards high yields. The overall increased performance accuracy in the low-yield regime additionally suggests that the balancing of the dataset improves the training of the model in general. However, it is important to note, that the 

## References:

[1] Rashmi Korlakai Vinayak, Ran Gilad-Bachrach. “DART: Dropouts meet Multiple Additive Regression Trees.” JMLR.[2] Yang, Yi, et al. "Multi-class active learning by uncertainty sampling with diversity maximization." International Journal of Computer Vision 113.2 (2015): 113-127.
[2] Yang, Yi, et al. "Multi-class active learning by uncertainty sampling with diversity maximization." International Journal of Computer Vision 113.2 (2015): 113-127.
[3] Branco, P., Torgo, L., Ribeiro, R. (2017). SMOGN: A Pre-Processing Approach for Imbalanced Regression. Proceedings of Machine Learning Research, 74:36-50. http://proceedings.mlr.press/v74/branco17a/branco17a.pdf.
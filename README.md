# Benchmark_Phospho_Identification
Deep learning-derived evaluation metrics enable effective benchmarking of computational tools for phosphopeptide identification

Tandem mass spectrometry (MS/MS)-based phosphoproteomics is a powerful technology for global phosphorylation analysis. However, applying four computational pipelines to a typical mass spectrometry (MS)-based phosphoproteomic dataset from a human cancer study, we observed a large discrepancy among the reported phosphopeptide identification and phosphosite localization results, underscoring a critical need for benchmarking. While efforts have been made to compare performance of computational pipelines using data from synthetic phosphopeptides, evaluations involving real application data have been largely limited to comparing the numbers of phosphopeptide identifications due to the lack of appropriate evaluation metrics. ![image](https://user-images.githubusercontent.com/48265327/130851731-75993b21-0536-4f7a-b651-e7341a283a8b.png)

We investigated three deep learning-derived features as potential evaluation metrics: phosphosite probability, Delta RT and spectral similarity. Predicted phosphosite probability is computed by MusiteDeep, which provides high accuracy as previously reported; Delta RT is defined as the absolute retention time (RT) difference between RTs observed and predicted by AutoRT; and spectral similarity is defined as the Pearson’s correlation coefficient between spectra observed and predicted by pDeep2. ![image](https://user-images.githubusercontent.com/48265327/130851777-df687692-4765-4081-a652-76da33ffe77d.png)

The benchmark metrics demonstrated in this study will enable users to select computational pipelines and parameters for routine analysis of phosphoproteomics data and will offer guidance for developers to improve computational methods.![image](https://user-images.githubusercontent.com/48265327/130851807-07cc61a4-add4-4d1f-ba38-0946b3239a1f.png)

![Figure1](https://user-images.githubusercontent.com/48265327/130851833-9fdf90b2-7a0e-434a-89c7-04a111cbec8a.png)

# Data Availability
The phosphoproteomics data that support the findings of this study are available in PRIDE (https://www.ebi.ac.uk/pride/)  with the identifier PXD01508746, PXD000138, PXD007145, PXD015284. Pre-prepared data for AutoRT base model training and testing were downloaded from the GitHub website (https://github.com/bzhanglab/AutoRT/tree/master/example/data). The phosphoproteomic data from uterine corpus endometrial carcinoma study (UCEC) were downloaded from the CPTAC data portal (https://cptac-data-portal.georgetown.edu/study-summary/S043). 
![image](https://user-images.githubusercontent.com/48265327/130851901-c88cc196-360d-48e1-ad22-fec94c4a6d4c.png)

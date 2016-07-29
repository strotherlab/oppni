# OPNI (Optimization of Preprocessing pipelines for NeuroImaging) for fMRI analysis by the Strother Laboratory

This repository provides OPNI (Optimization of Preprocessing pipelines for NeuroImaging), which does fast optimization of preprocessing pipelines for BOLD fMRI (Blood Oxygenation Level Dependent functional MRI). This software package identifies the set of preprocessing steps (“pipeline”) specific to each dataset, which optimizes quality metrics of Prediction and Reproducibility for post-processing analysis results (Strother et al., 2002;.This procedure has been shown to significantly improve signal detection, reliability of brain activations, and sensitivity to brain-behaviour correlations (Churchill et al., 2012a, 2012b). The pipeline software can also be used for simple automated batch-processing of fMRI datasets, if no appropriate analysis model is currently available to do optimization (e.g. some resting-state connectivity studies).

## What can OPNI do for you?

### Clean up your data: 
it will automatically output data that has been optimally processed to control for noise and artifact (e.g. due to motion, physiology, scanner noise), which you can then use for further analysis. Automated de-noising is done based on replicable, statistical criteria; no more tedious batch scripting, or hand-selection of ICA components!

### Analyze your data: 
as part of optimization, OPNI creates Z-scored maps of brain activity for each dataset. You can directly report these results, or take individual activation maps and do standard group-level analysis. OPNI can perform univariate or multivariate analysis of brain activity.  This can be done for a variety of different task designs, including event-related and block-design tasks. Recent additions also include seed-based connectivity and component modelling, i.e., selecting an optimal PCA subspace.

### Run batched pipelines: 
if you cannot analyze your data (or OPNI does not have the appropriate analysis tools), you can still process it using OPNI without doing optimization. Our scripts make it straightforward to choose the pipeline steps that you want, and perform automatic batch processing of large datasets.

## Implementation
Preprocessing and analysis scripts are coded in Matlab/Octave, and functions are called and managed using Python scripts. The code tests all possible combinations of a set of 12 different preprocessing steps, to identify the optimal pipeline for each fMRI dataset that is being tested. Current preprocessing options include AFNI utilities (Analysis of Functional NeuroImaging; Cox, 1996), along with a set of functions developed in-house. All steps are widely used in the fMRI literature, or demonstrated to be important in prior studies of pipeline optimization (e.g. Tegeler et al., 1999; La Conte et al., 2003; Shaw et al., 2003; Strother et al., 2004; Zhang et al., 2009; Churchill et al., 2012a, 2012b). Pipeline optimization is analysis-driven: it evaluates the quality of analysis results for each pipeline, via Prediction and Reproducibility metrics, and selects the pipeline that gives highest-quality outputs.  Analysis techniques are “modular” – you choose the task design and analysis model you wish to optimize, from a list of available models. The pipeline software also includes a procedure for automated spatial normalization of subjects to an anatomical template, using FSL utilities. This enables users to run group-level analysis of preprocessed results across subjects and task runs.

## References

Campbell K et al. (2013). Age Differences in the Intrinsic Functional Connectivity of Default Network Subsystems. Frontiers in Human Neuroscience. 5:73

Churchill NW et al. (2012a): Optimizing Preprocessing and Analysis Pipelines for Single-Subject fMRI: 2. Interactions with ICA, PCA, Task Contrast and Inter-Subject Heterogeneity. PLoS One. 7(2):e31147

Churchill NW et al. (2012b): Optimizing Preprocessing and Analysis Pipelines for Single-Subject FMRI. I. Standard Temporal Motion and Physiological Noise Correction Methods. Human Brain Mapping 33:609–627 

Churchill NW, Strother SC (2013): PHYCAA+: an optimized, adaptive procedure for measuring and controlling physiological noise in BOLD fMRI. NeuroImage 82:306-325.

Cox RW (1996): AFNI: Software for analysis and visualization of functional magnetic resonance neuroimages. Computers and Biomedical Research, an International Journal, 29(3): 162-173. 

Glover GH, et al. (2000): Image-based method for retrospective correction of physiological motion effects in fMRI: RETROICOR. Magnetic Resonance in Medicine: Official Journal of the Society of Magnetic Resonance in Medicine / Society of Magnetic Resonance in Medicine. 44(1): 162-167. 

LaConte S et al. (2003): The Evaluation of Preprocessing Choices in Single-Subject BOLD fMRI Using NPAIRS Performance Metrics. NeuroImage, 18(1):10-27

Schmah T, et al. (2010): Comparing classification methods for longitudinal fMRI studies. Neural Comput, 22:2729-62.

Shaw ME et al. (2003): Evaluating subject specific preprocessing choices in multisubject fMRI data sets using data-driven performance metrics. NeuroImage.  19(3):988-1001.

Strother SC, Anderson J, Hansen LK, Kjems U, Kustra R et al. (2002): The Quantitative Evaluation of Functional Neuroimaging Experiments: The NPAIRS Data Analysis Framework. NeuroImage 15:747–771

Strother S et al. (2004): Optimizing the fMRI data-processing pipeline using prediction and reproducibility performance metrics: I. A preliminary group analysis. NeuroImage. 23:S196-S207.

Strother S et al. (2014): Stability and Reproducibility in fMRI Analysis, in Practical Applications of Sparse Modeling, I. Rish, G. A. Cecchi, A. Lozano, and A. Niculescu-Mizil, Eds., pp. 99-121, Boston: MIT Press.

Tegeler C et al. (1999): Reproducibility of BOLD-based functional MRI obtained at 4 T. Hum. Brain Mapp. 7:267–283

Zhang J et al. (2009). Evaluation and optimization of fMRI single-subject processing pipelines with NPAIRS and second-level CVA. Magn. Reson. Imag. 27:264–278


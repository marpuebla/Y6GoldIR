# Y6GoldIR
Enhancing the Dark Energy Survey Y6 Gold catalog redshifts using infrared data.

This repository contains the code and analysis developed to estimate **photometric redshifths (photo-z)** for galaxies from the **Dark Energy Survey (DES)** by combining optical data with infrared measurements from 
the **Wide-field Infrared Survey Explorer (WISE)**.

This work is based on the **Angular Neighbourhood Fitting (ANF)** machine learning algorithm, eveloped by De Vicente et al. (2016), and evaluates how the inclusion of infrared photometry affects photo-z performance in terms of **bias**, **scatter** and **outlier fraction**. 

## ANF:

This software, developed in Python, allows for the calculation of photometric redshifts by applying advanced fitting techniques. In ANF, the photometric redshift of a galaxy is calculated based on the redshifts of nearby training galaxies within the multimagnitude space. The redshifts of the training sample define a hypersurface across the multimagnitude space. In this context, for any given point in this space, which corresponds to a galaxy with unknown redshift, the shape of the local hypersurface can vary in smoothness depending on the direction. Accurate redshift estimation is more likely when the neighbourhood is selected along a direction in which the redshift hypersurface varies smoothly, allowing for a better fit to the local structure. 

ANF algorithm provide the main photo-z value as: DNF_Z, determined by the fit of a number of neighbour galaxies to a hyperplane in the magnitude space. The obtained photometric redshift for each object is called (DNF_Z), determined using around 80 neighbours.

## The metrics:

The sample we are analysing includes magnitude measurements from the broad filters of all three surveys, as well as the spectroscopic redshifts. Taking the spectroscopic measurement as the truth value for the redshift, we can evaluate the quality of the photometric redshift estimation. 

1.Bias, $\mu$: measures the systematic offset between $z_{phot}$ and $z_{spec}$, where $z_{phot}$ and $z_{spec}$ refer to the photometric and spectroscopic redshift, respectively.
The individual bias for each galaxy is defined as:

Δ
𝑧
=
𝑧
𝑝
ℎ
𝑜
𝑡
−
𝑧
𝑠
𝑝
𝑒
𝑐
Δz=z
phot
	​

−z
spec
	​


If we consider the number of photometric redshift galaxies as $N$, and the mean of the absolute values, the metric is set as:

𝜇
=
1
𝑁
∑
𝑖
=
1
𝑁
∣
Δ
𝑧
𝑖
∣
μ=
N
1
	​

i=1
∑
N
	​

∣Δz
i
	​

∣

A small value of this metric indicates that photometric redshifts are accurate with respect to spectroscopic redshifts.

We used the bootstrap  resampling statistical method to estimate the uncertainties and error bars for the metrics. Bootstrap is a non parametric technique that allows the estimation of the standard error and the confidence intervals of a given statistic by repeatedly resampling the original data with replacement. We implemented the bootstrap method using the SciPy library from Python. For each metric (the mean bias, scatter or outlier fraction), we generated 100 bootstrap resamples for the original matched catalogue. The metric is recomputed on each resample, and the resulting distribution was used to estimate the standard error. We adopted a 95$\%$ confidence level, and the resulting standard error is used as the symmetric uncertainty for the corresponding bin in our analysis.

---
## Examples 

Metrics for the match between DES and WISE: 

a) Photometric vs spectroscopic redshift

<img width="400" height="500" alt="F_1_b_wise" src="https://github.com/user-attachments/assets/ec2da19b-3b88-4567-86a1-094af47035eb" />

b) Bias

<img width="400" height="500" alt="F_8_b_wise" src="https://github.com/user-attachments/assets/d1d98bb9-b0df-4cbf-a693-3a5a4d671966" />

c) The 68% quantile error over (1 + z)

<img width="400" height="500" alt="F_4_b_wise" src="https://github.com/user-attachments/assets/b4ca11f2-ba90-4aa3-870a-50dec8e5fecf" />

d) Outliers

<img width="400" height="500" alt="F_7_banerji_b_wise" src="https://github.com/user-attachments/assets/05842e8f-a82a-46f7-b57d-28c29b6c6b91" />

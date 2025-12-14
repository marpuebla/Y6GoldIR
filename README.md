# Y6GoldIR
Enhancing the Dark Energy Survey Y6 Gold catalog redshifts using infrared data.

This repository contains the code and analysis developed to estimate **photometric redshifths (photo-z)** for galaxies from the **Dark Energy Survey (DES)** by combining optical data with infrared measurements from 
the **Wide-field Infrared Survey Explorer (WISE)**.

This work is based on the **Angular Neighbourhood Fitting (ANF)** machine learning algorithm, developed by De Vicente et al. (2016), and evaluates how the inclusion of infrared photometry affects photo-z performance in terms of **bias**, **scatter** and **outlier fraction**. 

## ANF:

This software, developed in Python, allows for the calculation of photometric redshifts by applying advanced fitting techniques. In ANF, the photometric redshift of a galaxy is calculated based on the redshifts of nearby training galaxies within the multimagnitude space. The redshifts of the training sample define a hypersurface across the multimagnitude space. In this context, for any given point in this space, which corresponds to a galaxy with unknown redshift, the shape of the local hypersurface can vary in smoothness depending on the direction. Accurate redshift estimation is more likely when the neighbourhood is selected along a direction in which the redshift hypersurface varies smoothly, allowing for a better fit to the local structure. 

ANF algorithm provide the main photo-z value as: DNF_Z, determined by the fit of a number of neighbour galaxies to a hyperplane in the magnitude space. The obtained photometric redshift for each object is called (DNF_Z), determined using around 80 neighbours.

## The metrics:

The sample we are analysing includes magnitude measurements from the broad filters of all three surveys, as well as the spectroscopic redshifts. Taking the spectroscopic measurement as the truth value for the redshift, we can evaluate the quality of the photometric redshift estimation. 

We used the bootstrap  resampling statistical method to estimate the uncertainties and error bars for the metrics. Bootstrap is a non parametric technique that allows the estimation of the standard error and the confidence intervals of a given statistic by repeatedly resampling the original data with replacement. We implemented the bootstrap method using the SciPy library from Python. For each metric (the mean bias, scatter or outlier fraction), we generated 100 bootstrap resamples for the original matched catalogue. The metric is recomputed on each resample, and the resulting distribution was used to estimate the standard error. We adopted a 95$\%$ confidence level, and the resulting standard error is used as the symmetric uncertainty for the corresponding bin in our analysis.


**Bias,** $\mu$: measures the systematic offset between $z_{phot}$ and $z_{spec}$, where $z_{phot}$ and $z_{spec}$ refer to the photometric and spectroscopic redshift, respectively.  
The individual bias for each galaxy is defined as:

$$
\Delta z = z_{phot} - z_{spec}
$$

If we consider the number of photometric redshift galaxies as $N$, and the mean of the absolute values, the metric is set as:

$$
\mu = \frac{1}{N}\sum_{i=1}^{N} |\Delta z_i|
$$

A small value of this metric indicates that photometric redshifts are accurate with respect to the spectroscopic redshifts.

---

**Dispersion,** $\sigma$: represents the scatter of the individual photometric residuals ($\Delta z_i$) for each galaxy around their mean value.  
It quantifies the variability of the photometric redshift errors across the sample.

If we consider:

$$
b = \frac{1}{N}\sum_{i=1}^{N} (\Delta z_i)
$$

then the dispersion is set as:

$$
\sigma = \sqrt{\frac{1}{N}\sum_{i=1}^{N} (\Delta z_i - b)^2}
$$

---

**Precision in 68-quantile,** $\sigma_{68}$: quantifies the scatter or dispersion of the photometric redshifts.  
It represents the width of the distribution of photometric redshifts around the median that contains 68% of the data points.  
Specifically, it corresponds to the 68% quantile error, and is defined as:

$$
\sigma_{68} = \frac{1}{2}(P_{84} - P_{16})
$$

where $P_{84}$ and $P_{16}$ are the 84th and 16th percentiles of the cumulative distribution, respectively.  
The value of $\sigma_{68}$ is indicative of how well the photometric redshift estimates cluster around the spectroscopic redshift.  
A small value suggests that the estimation is more accurate.

---

**Normalised** $\sigma_{68}$: the 68% quantile error over $(1 + z)$, computed as:

$$
\frac{\sigma_{68}}{1 + z_{spec}} = \frac{1}{2}\frac{P_{84} - P_{16}}{1 + z_{spec}}
$$

The normalisation is particularly useful for ensuring that the error does not systematically increase with redshift; otherwise, it would be a photo-$z$ estimation bias.

---

**Outlier fraction:** this metric quantifies the fraction of objects with large biases, for which the photometric redshift estimates are significantly inaccurate.  
It is defined following the criterion of *Banerji et al. (2015)*, where an object is classified as an outlier if the following condition is satisfied:

$$
\frac{\Delta z}{1 + z_{spec}} > 0.15
$$

This fixed-threshold definition is preferred over $\sigma_{68}$-based criteria ($\Delta z / (1 + z_{spec}) > n \cdot \sigma_{68}$), as the value of $\sigma_{68}$ tends to increase with redshift.  
Therefore, applying a threshold that scales with $\sigma_{68}$ can lead to an underestimation of the outlier rate in higher redshift bins.  
In contrast, the Banerji criterion sets a constant threshold for all redshift bins, allowing for stronger outlier identification and a direct comparison across the entire redshift range.




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

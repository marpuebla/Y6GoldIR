# Y6GoldIR
Enhancing the Dark Energy Survey Y6 Gold catalog redshifts using infrared data.

This repository contains the code and analysis developed to estimate **photometric redshifths (photo-z)** for galaxies from the **Dark Energy Survey (DES)** by combining optical data with infrared measurements from 
the **Wide-field Infrared Survey Explorer (WISE)**.

This work is based on the **Angular Neighbourhood Fitting (ANF)** machine learning algorithm, and evaluates how the inclusion of infrared photometry affects photo-z performance in terms of **bias**, **scatter** and **outlier fraction**. 

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

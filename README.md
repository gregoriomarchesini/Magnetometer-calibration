# Magnetometer calibration techniques

# Introduction
The following repository contains three different scripts for the calibration of a tri-axes magnetometer which needs to be calibrated for hard/soft ferromagnetic interfernce and scale factors .
The three calibration procedures proposed herein are **Manual Calibration(MC)**,**Non-Linear Least Square Calibration (NLLS)** and **Adjusted Ordinary Least Square Calibration (ALQ)**.

# Summary 
Here is brief summary of the functions that can be found in the repository for completeness 

|LIST OF FUNCTIONS|Brief Description|
|---------------- |---|
|```ALS_calibration.m.```        | ALS Calibration script|
|```Calibration_data.mat. ```    | File where you can same your solved parameters after each calibartion|
|```First_order_calibration.m``` | MC calibration script|
|```fitter.m        ```          | NLLS calibration script|
|```mag_batch.mat  ```           | Batch of random measurments to use as samples for the calibration algorithms|
|```vec_s_inv.m   ```            | Auxiliary routine needed for ALS calibration|
|``` vec_s.m     ```             | Auxiliary routine needed for ALS calibration|


# References

[1] Markovsky, Ivan & Kukush, Alexander & Huffel, Sabine. (2004).[*Consistent least squares fitting of ellipsoids*](https://www.researchgate.net/publication/39994614_Consistent_least_squares_fitting_of_ellipsoids). Numerische Mathematik. 98. 10.1007/s00211-004-0526-9. 

[2] Renaudin, Valérie, Muhammad Haris Afzal, and Gérard Lachapelle. [*Complete triaxis magnetometer calibration in the magnetic domain.*](https://www.researchgate.net/publication/303721929_Complete_triaxis_magnetometer_calibration_in_the_magnetic_domain) Journal of sensors 2010 (2010).

[3] Kelley, Carl T. Iterative methods for optimization. Society for Industrial and Applied Mathematics, 1999.


May the Force be with you 

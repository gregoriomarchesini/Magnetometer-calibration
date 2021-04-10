# Magnetometer_calibration_techniques

# Introduction
The following repository contains three different script for the calibration of a tri-axes magnetometer which needs to be calibrated for hard/soft ferromagnetic interfernce and scale factors .
The three calibration procedure proposed herein are **Manual Calibration(MC)**,**Non-Linear Least Square Calibration (NLLS)** and **Adjusted Ordinary Least Square Calibration (ALQ)**.
All the three methods are briefly explained here together with a complete description of how to use the given scripts to solve the calibration for each method.

It is important to note that only scale factor and hard ferromagnetic interfernce will be considered in the calibration as main sources of measurement bias, although non-orthogonalitoes 
and soft ferromagnetic interference are present. The only method that is able to calibrate all the aforementioned sources of bias (hard/soft ferromagnetic interfernce, scale factors and 
non-orthogonalities) is the ALQ method, which development is based on [3]. Both MC and NLLS solve only for hard ferromagnetic inetrfernce (offset from now on) and scale factor.

# Manual Calibration

# LLSQ Calibration 

# Adjuasted Least Square Calibration






# References

[1] Markovsky, Ivan & Kukush, Alexander & Huffel, Sabine. (2004).* Consistent least squares fitting of ellipsoids*. Numerische Mathematik. 98. 10.1007/s00211-004-0526-9. 

[2] Renaudin, Valérie, Muhammad Haris Afzal, and Gérard Lachapelle. *Complete triaxis magnetometer calibration in the magnetic domain.* Journal of sensors 2010 (2010).

[3] Kelley, Carl T. Iterative methods for optimization. Society for Industrial and Applied Mathematics, 1999.

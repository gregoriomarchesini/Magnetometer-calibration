# Magnetometer calibration techniques

# Introduction
The following repository contains three different script for the calibration of a tri-axes magnetometer which needs to be calibrated for hard/soft ferromagnetic interfernce and scale factors .
The three calibration procedure proposed herein are **Manual Calibration(MC)**,**Non-Linear Least Square Calibration (NLLS)** and **Adjusted Ordinary Least Square Calibration (ALQ)**.
All the three methods are briefly explained here together with a complete description of how to use the given scripts to solve the calibration for each method.

It is important to note that only scale factor and hard ferromagnetic interfernce will be considered in the calibration as main sources of measurement bias, although non-orthogonalitoes 
and soft ferromagnetic interference are present. The only method that is able to calibrate all the aforementioned sources of bias (hard/soft ferromagnetic interfernce, scale factors and 
non-orthogonalities) is the ALQ method, which development is based on [3]. Both MC and NLLS solve only for hard ferromagnetic inetrfernce (offset from now on) and scale factor.

Given the this semplifications, it is important to define a simple model that describes how the magnetometer senses and external magnetic field. A complete and precise mathematical description of a magnetometer can be found in [3], but a simplified version is given herein for pratical reasons. Each axis of the magnetometer gives one of the three componet of a refernce magnetic field that it is measured (it can be an artificial magnetic field or the geomagnetic fiedl in most of the cases). However, the measurnments along each axis is biased by a ceratain offset (hard ferromagnetic interference) and a ceratin scale factor. This means in the former case that even with zero external field, the magnetomegter senses a nn-zero magnetic field that is diffrent for each axis. In the latter case, the problem is that the external magnetic field is measured with a constant multiplicative factor that is different for each axis. Hence it is possible to describe the measured magnetic field along each axis (x axis of the magnetometer for example) as:

H<sup>sensed</sup> = (H<sup>real</sup> + **off**) x **s**
  
Where H<sup>real</sup> is the true magnetic field componet along one axis, **off** is the offset and **s** is the scale factor. The whole iam of the calibration is to find the value of offset and scale factor for each axis, so that at the end it will be possible to obtain correct measurements from the magnetometer using the obtained calibration parameters. A total of six parameters are to be determined that can be stacked into a vector 

beta = [off<sub>x</sub> off<sub>y</sub> off<sub>z</sub> s<sub>x</sub> s<sub>y</sub> s<sub>z</sub>]

Graphically this correspond in resetting the measurements from an ellipsoid to a sphere with radious equal to the refernce magnetic fiedl intensity. This is becasuse, the magnetic field must always have the same intensity even if the magnetometer is rotating if no bias is affecting the measurements or if the bias is properly calibrated. The compoent of the magnetic field will change as the magnetometer changes orientation, but ideally the magnitude of the magnetic field vecotor must be constant.

![Skectch](images/Mist_scatter.png)
  
# Manual Calibration

## Description
 
The MC coalibration consist in two simple steps that are to be repeated for each axis, for a total of three times. First it will be necessary to have at disposal a known magnetic field field in both direction and intesity. This can be an external artificial magnetic field or or the geomagnetic field which can be obtained from <code>igrf()</code> in Matlab for example (or an online calculator like the one [here](http://www.geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html).

Once the refernce is derfined the only task to accomplish is pointing the magnetometer axis aligned with the refernce magnetic field. Then one measurements should be taken in the same direction of the reference and a second one should be taken in the opposit direction. The following figure illustartes the process :


![Skectch](images/new_cal.png)

From this two measurements it is possible to solve for **off** and **s** at each axis axis using this symple system of equations 

H<sup>sensed(+)</sup> = (H<sup>real</sup> + **off**) x **s**
H<sup>sensed(-)</sup> = (-H<sup>real</sup> + **off**) x **s**

where H<sup>real</sup> is known and H<sup>sensed(-)</sup> and H<sup>sensed(+)</sup> are derived from the measurments (in the same and opposit direction of the reference respectively)








# LLSQ Calibration 

# Adjuasted Least Square Calibration






# References

[1] Markovsky, Ivan & Kukush, Alexander & Huffel, Sabine. (2004).* Consistent least squares fitting of ellipsoids*. Numerische Mathematik. 98. 10.1007/s00211-004-0526-9. 

[2] Renaudin, Valérie, Muhammad Haris Afzal, and Gérard Lachapelle. *Complete triaxis magnetometer calibration in the magnetic domain.* Journal of sensors 2010 (2010).

[3] Kelley, Carl T. Iterative methods for optimization. Society for Industrial and Applied Mathematics, 1999.

CC-CSI
======

Introduction
------------

CC-CSI is a MATLAB-based package of the **Cross-Correlated Contrast Source Inversion** algorithm. This package inverts the [Fresnel data](http://www.fresnel.fr/3Ddatabase/) with multiple frequency components. 

- Basic Formulation

	- data error equation, 

	![](http://latex.codecogs.com/gif.latex?%5Crho_%7Bp%2Ci%7D%20%3D%20y_%7Bp%2Ci%7D%20-%20%5CPhi_%7Bp%2Ci%7Dj_%7Bp%2Ci%7D)
	$$ $$

	- state error equation, 

	![](http://latex.codecogs.com/gif.latex?%5Cgamma_%7Bp%2Ci%7D%20%3D%20%5Cchi_i%20e%5E%7B%5Ctext%7Binc%7D%7D_%7Bp%2Ci%7D%20&plus;%20%5Cchi_i%20A_i%5E%7B-1%7Dj_%7Bp%2Ci%7D%20-%20j_%7Bp%2Ci%7D)

	- cross-correlated error equation, 

	![](http://latex.codecogs.com/gif.latex?%5Cxi_%7Bp%2Ci%7D%20%3D%20y_%7Bp%2Ci%7D%20-%20%5CPhi%5Cleft%28%5Cchi_i%20e%5E%7B%5Ctext%7Binc%7D%7D_%7Bp%2Ci%7D%20&plus;%20%5Cchi_i%20A_i%5E%7B-1%7Dj_%7Bp%2Ci%7D%5Cright%29)

	- CC-CSI cost functional, 

	![](http://latex.codecogs.com/gif.latex?%5Csum_i%5Cleft%28%20%5Ceta_i%5E%5Cmathcal%7BS%7D%5Csum_p%5Cleft%5C%7C%5Crho_%7Bp%2Ci%7D%5Cright%5C%7C%5E2_%5Cmathcal%7BS%7D%20&plus;%20%5Ceta_i%5E%5Cmathcal%7BD%7D%5Csum_p%5Cleft%5C%7C%5Cxi_%7Bp%2Ci%7D%5Cright%5C%7C%5E2_%5Cmathcal%7BD%7D%20&plus;%20%5Ceta_i%5E%5Cmathcal%7BS%7D%5Csum_p%5Cleft%5C%7C%5Cxi_%7Bp%2Ci%7D%5Cright%5C%7C%5E2_%5Cmathcal%7BS%7D%5Cright%29)

	$$\sum_i\left( \eta_i^\mathcal{S}\sum_p\left\|\rho_{p,i}\right\|^2_\mathcal{S} + \eta_i^\mathcal{D}\sum_p\left\|\xi_{p,i}\right\|^2_\mathcal{D} + \eta_i^\mathcal{S}\sum_p\left\|\xi_{p,i}\right\|^2_\mathcal{S}\right)$$

This package also contains the MR-CSI (multiplicative-regularized CSI) and CSI algorithms for comparison.

- See `INSTALL.md` for installation instruction.

- For details of CC-CSI, see

	- [*Cross-Correlated Contrast Source Inversion*](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7862846), S. Sun, B. J. Kooij, T. Jin, and A. G. Yarovoy, IEEE Transactions on Antennas and Propagation, 65 (5), 2592 - 2603, 2017

	- [*Inversion of Multifrequency Data With the Cross‐Correlated Contrast Source Inversion Method*](https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2017RS006505), S. Sun, B. J. Kooij, and A. G. Yarovoy, Radio Science, 53 (6), 710-723, 2018


Feedback
--------
We would be delighted to hear from you if you find CC-CSI useful, or if you have any suggestions, contributions, or bug reports. Please send these to

Shilong Sun (shilongsun@icloud.com)




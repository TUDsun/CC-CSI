CC-CSI
======

Introduction
------------

CC-CSI is a MATLAB-based package of the **Cross-Correlated Contrast Source Inversion** algorithm. This package inverts the [Fresnel data](http://www.fresnel.fr/3Ddatabase/) with multiple frequency components. 

- The error terms 

	- Data error

	<!--- http://latex.codecogs.com/eqneditor/editor.php -->

	<!--- \xi_{p,i} = y_{p,i} - \Phi_{p,i}\left(\chi_i e^{\text{inc}}_{p,i} + \chi_i A^{-1}_ij_{p,i}\right) -->

	<div align=center><img src="http://latex.codecogs.com/gif.latex?%5Crho_%7Bp%2Ci%7D%20%3D%20y_%7Bp%2Ci%7D%20-%20%5CPhi_%7Bp%2Ci%7Dj_%7Bp%2Ci%7D"/></div>

	- State error

	<!--- \gamma_{p,i} = \chi_i e^{\text{inc}}_{p,i} + \chi_i A_i^{-1}j_{p,i} - j_{p,i}  -->

	<div align=center><img src="http://latex.codecogs.com/gif.latex?%5Cgamma_%7Bp%2Ci%7D%20%3D%20%5Cchi_i%20e%5E%7B%5Ctext%7Binc%7D%7D_%7Bp%2Ci%7D%20&plus;%20%5Cchi_i%20A_i%5E%7B-1%7Dj_%7Bp%2Ci%7D%20-%20j_%7Bp%2Ci%7D"/></div>


	- Cross-correlated error 

	<!--- \xi_{p,i} = y_{p,i} - \Phi_{p,i}\left(\chi_i e^{\text{inc}}_{p,i} + \chi_i A^{-1}_ij_{p,i}\right) -->

	<div align=center><img src="http://latex.codecogs.com/gif.latex?%5Cxi_%7Bp%2Ci%7D%20%3D%20y_%7Bp%2Ci%7D%20-%20%5CPhi_%7Bp%2Ci%7D%5Cleft%28%5Cchi_i%20e%5E%7B%5Ctext%7Binc%7D%7D_%7Bp%2Ci%7D%20&plus;%20%5Cchi_i%20A%5E%7B-1%7D_ij_%7Bp%2Ci%7D%20%5Cright%20%29"/></div>

- The cost functional of multi-frequency CC-CSI

	<!--- \mathcal{C}_{\text{CC-CSI}} = \sum_i\eta^\mathcal{S}_i\sum_p\left\|\rho_{p,i}\right\|^2_\mathcal{S}+\sum_i\eta^\mathcal{D}_i\sum_p\left\|\gamma_{p,i}\right\|^2_\mathcal{D}+\sum_i\eta^\mathcal{S}_i\sum_p\left\|\xi_{p,i}\right\|^2_\mathcal{S} -->

	<div align=center><img src="http://latex.codecogs.com/gif.latex?%5Cmathcal%7BC%7D_%7B%5Ctext%7BCC-CSI%7D%7D%20%3D%20%5Csum_i%5Ceta%5E%5Cmathcal%7BS%7D_i%5Csum_p%5Cleft%5C%7C%5Crho_%7Bp%2Ci%7D%5Cright%5C%7C%5E2_%5Cmathcal%7BS%7D&plus;%5Csum_i%5Ceta%5E%5Cmathcal%7BD%7D_i%5Csum_p%5Cleft%5C%7C%5Cgamma_%7Bp%2Ci%7D%5Cright%5C%7C%5E2_%5Cmathcal%7BD%7D&plus;%5Csum_i%5Ceta%5E%5Cmathcal%7BS%7D_i%5Csum_p%5Cleft%5C%7C%5Cxi_%7Bp%2Ci%7D%5Cright%5C%7C%5E2_%5Cmathcal%7BS%7D"/></div>

This package also contains the MR-CSI (multiplicative-regularized CSI) and CSI algorithms for comparison.

- See `INSTALL.md` for installation instruction.


- We would be delighted to hear from you if you find CC-CSI useful, or if you have any suggestions, contributions, or bug reports. Please send these to 

	- Shilong Sun (shilongsun@icloud.com)


Citation
--------

- [*Cross-Correlated Contrast Source Inversion*](https://arxiv.org/abs/1906.10864), S. Sun, B. J. Kooij, T. Jin, and A. G. Yarovoy, IEEE Transactions on Antennas and Propagation, 65 (5), 2592 - 2603, 2017

- [*Inversion of Multifrequency Data With the Cross‐Correlated Contrast Source Inversion Method*](https://arxiv.org/abs/1906.10814), S. Sun, B. J. Kooij, and A. G. Yarovoy, Radio Science, 53 (6), 710-723, 2018




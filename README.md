# ComplexBessel
 Complex Bessel Function Numerical Evaluation Note
 
## Im(z) &lt; 0
Complex Bessel functions are analytically connected so that they are complex conjugate symmetrically about the real axis.

![bessel conj](figures/bessel_conj.svg)  

[DLMF 10.11](https://dlmf.nist.gov/10.11)  
[DLMF 10.34](https://dlmf.nist.gov/10.34)  
 
## Re(z) &lt; 0, Im(z) &geq; 0
Complex Bessel functions are defined to be discontinuous on the negative real axis.

![besselj minus rez](figures/besselj_minus_rez.svg)  
![bessely minus rez](figures/bessely_minus_rez.svg)  
![besseli minus rez](figures/besseli_minus_rez.svg)  
![besselk minus rez](figures/besselk_minus_rez.svg)  
![hankel1 minus rez](figures/hankel1_minus_rez.svg)  
![hankel2 minus rez](figures/hankel2_minus_rez.svg)  

[DLMF 10.11](https://dlmf.nist.gov/10.11)  
[DLMF 10.34](https://dlmf.nist.gov/10.34)  

## Re(z) &geq; 0, Im(z) &geq; 0, Asymptotic Expansion

![hankel coef](figures/hankel_coef.svg)  
![hankel omega](figures/hankel_omega.svg)  

![besselj asymp](figures/besselj_asymp.svg)  
![bessely asymp](figures/bessely_asymp.svg)  
![besseli asymp](figures/besseli_asymp.svg)  
![besselk asymp](figures/besselk_asymp.svg)  
![hankel1 asymp](figures/hankel1_asymp.svg)  
![hankel2 asymp](figures/hankel2_asymp.svg)  

[DLMF 10.17](https://dlmf.nist.gov/10.17)  
[DLMF 10.40](https://dlmf.nist.gov/10.40)  

## Numerical Evaluation Stability

|BesselJ, BesselY|BesselI|BesselK|
|---|---|---|
|![besseljy convergence](figures/besseljy_convergence.svg)|![besseli convergence](figures/besseli_convergence.svg)|![besselk convergence](figures/besselk_convergence.svg)|

## Connection Formula

![besselji](figures/besselji.svg)  
![besselyk](figures/besselyk.svg)  
![bessel itoj](figures/bessel_itoj.svg)  

[DLMF 10.27](https://dlmf.nist.gov/10.27)  
 
## Reference
[DLMF](https://dlmf.nist.gov/10)  
[Wolfram Math World](https://mathworld.wolfram.com/BesselFunction.html)

## Author

[T.Yoshimura](https://github.com/tk-yoshimura)

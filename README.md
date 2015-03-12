Full Euler solver on deep water
=======================

The present script is a Fourier-type pseudo-spectral solver of the full Euler equations with free surface on a fluid layer of infinite depth (the so-called *deep water approximation*). The time-dependent fluid domain is transformed into a strip using the conformal mapping technique (coined by L.V. Ovsyannikov for the first time and developed later by A. Dyachenko *et al.*).

## Acknowledgements

The author would like to thank the following colleagues (in alphabetical order) who helped him to develop and initialize correctly the present solver:

* [Didier Clamond](http://math.unice.fr/~didierc/), [University of Nice Sophia Antipolis](http://unice.fr/), [Laboratoire J.A. Dieudonné](http://math.unice.fr/) France

* [Bernard Ee](https://www.researchgate.net/profile/Bernard_Ee2), [Tel Aviv University](http://english.tau.ac.il/), [School of Mechanical Engineering](http://engineering.tau.ac.il/), Israel

* [Lev Shemer](http://www.eng.tau.ac.il/~shemer/), [Tel Aviv University](http://english.tau.ac.il/), [School of Mechanical Engineering](http://engineering.tau.ac.il/), Israel
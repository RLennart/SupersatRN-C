# Supersaturation Rate Numerical-Calculator
<p align="center">
<picture>
 <source media="(prefers-color-scheme: dark)" srcset="logo_supersatrnc.svg">
 <source media="(prefers-color-scheme: light)" srcset="logo_supersatrnc.svg">
 <img alt="Logo SupersatRN-C" src="logo_supersatrnc.svg">
</picture>
</p>
Supersaturation Rate Numerical Calculator (SupersatRN-C) is a simple software package for simulating mass transfer in thin-film solution printing of perovskite absorber layers. It allows the calculation of thickness transients of solution films for different methods of increasing the supersaturation, so-called quenching techniques. The core is an ODE solver working on multicomponent mass trasnfer equations introduced in the publication [ref]. The mass transfer is determined Sherwood correlations from textbooks - distilled in 'correlations.py'. The material properties of the perovskite solutions and their respective literature sources are encoded in the file 'material_properties.py' and defined as classes in 'materials.py'. The software package is provided with the funcionality of the two solution formulations Methylammonium Lead Iodide (MAPI) dissolved in N,N-Dimethylformamide and Double Cation Perovsktie (DCP) dissolved in Dimethylsulfoxide (DMSO):DMF:Gamma-Butyrolactone (GBL). However, it is deliberately designed to be extensible and modular, so that more solution formulations can be added in the future. The principal output parameter of the SupersatRN-C is the supersaturation rate calculated by the temporal derivative of 
ln(C/C0)  at the point in time where the film reaches critical concentration. 

### Basic Features
The basic features of SupersatRN-C are 

* Calculation of dynamic evolution of film thickness part of every component in a perovskite solution film for different solutions and different processing methods.
* Systematic treatment of necessary processing parameters for full specification, providing a transparent documentation of their impact on the mass transfer dynamics.
* Determination of supersaturation rates of all common quenching methods used in perovskite solution processing, enabling quantitative comparison between these.
* Intuitive visualization of mass transfer dynamics along with a plot of supersaturation, supersaturation rate and optionally interferometric data.
* Possibility of regime-based processing methods combining arbitrary methods and/or activity changes within the solution.
* Possibility to fit activity coefficients based on interferometric data.

### Credits 
This software project was intiated and impelmented by **[Simon Ternes](mailto:ternes@ing.uniroma2.it?subject=[GitHub]%20Question%20on%SupersatRN-C)** under the supervision of **[Ulrich W. Paetzold](mailto:ulrich.paetzold@kit.edu)**.

The financial support by the following **projects and grants** is gratefully acknowledged:

- [PERCISTAND](https://percistand.eu/en) (funding code: 850937), European Union's Horizon 2020 research and innovation programme
- Helmholtz Young Investigator Group of U. W. Paetzold (funding code: VH-NG-1148), [Helmholtz Association](https://www.helmholtz.de/)
- [PEROSEED](https://www.helmholtz-berlin.de/projects/peroseed/index_en.html) (funding code: ZT-0024), [Helmholtz Association](https://www.helmholtz.de/)
- CAPITANO (funding code: 03EE1038B), [Federal Ministry for Economic Affairs and Energy](https://www.bmwi.de/)
- 27Plus6 (funding code: 03EE1056B), [Federal Ministry for Economic Affairs and Energy](https://www.bmwi.de/)

### Contributing

If you want to contribute to this project and make it better, your help is very welcome! 

### Contact

For any questions regarding the software, please contact **[Simon Ternes](mailto:ternes@ing.uniroma2.it?subject=[GitHub]%20Question%20on%SupersatRN-C)**. See also [here](http://www.chose.uniroma2.it/en/postdoc/421-simon-ternes.html).

### Citing

If you use our software or parts of it in the current or a modified version, you are obliged to provide proper attribution. This can be to our paper describing the software:
[ref]

### License

This software is licensed under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html) license. © 2021 EYcalc - Simon Ternes

Interested in a sublicense agreement to use SupersatRN-C in a non-free/restrictive environment? Contact **[Simon Ternes](mailto:ternes@ing.uniroma2.it?subject=[GitHub]%20Question%20on%SupersatRN-C)**!



# ANSYS Fatigue Cohesive Element
This code repository contains a user element subroutine for an ANSYS Mechanical APDL cohesive element to model 2D delamination onset/initiation and propagation under quasi-static and fatigue loading. The user element calculates also the J-integral in the cohesive zone and can be used in models with multiple delaminations. Tests cases showcasing these capabilities are also included in this repository.

The code was developed at the CraCS research group at Aalborg University, Siemens Gamesa RE and University of Girona by Iñigo Urcelay Oca, Brian L.V. Bak, Albert Turon, and Esben Lindgaard.

LICENSE
---------------
Copyright (c) 2023, CraCS Research Group at Aalborg University, Siemens Gamesa RE, and Universitat de Girona 
All rights reserved.

This project is licensed under the terms of the GNU Lesser General Public License v3.0 license.

The author would like to kindly emphasize that the following rules are respected when using the code. 
* Use of the code in source and binary forms, with or without modification are permitted provided that due reference is given to the author, the code and publications describing the implemented models.
* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

CONTENT
---------------
This code repository contains the following files:
* USERELEM_library_files: ANSYS files created after compiling and linking a user-element subroutine through the method of creating a dynamic library
* Static_DD_DCB.inp: ANSYS Mechanical APDL input file to model a Double-Delamination Double Cantilever Beam (DD DCB) specimen under quasi-static loading
* Fatigue_DCB.inp: ANSYS Mechanical APDL input file to model a Double Cantilever Beam (DCB) specimen under fatigue loading
* Fatigue_ENF.inp: ANSYS Mechanical APDL input file to model a End Notched Flexure (ENF) specimen under fatigue loading
* Fatigue_MMB.inp: ANSYS Mechanical APDL input file to model a Mixed Mode Bending (MMB) specimen under fatigue loading
* Fatigue_DD_DCB.inp: ANSYS Mechanical APDL input file to model a Double-Delamination Double Cantilever Beam (DD DCB) specimen under fatigue loading
* Fatigue_DNS.inp: ANSYS Mechanical APDL input file to model a Double-Notched Shear (DNS) specimen under fatigue loading

HOW TO USE
---------------
Any ANSYS subroutine source code must be first compiled and linked to ANSYS. There exist three different methods to do so. Refer to ANSYS own documentation for a more detailed explanation. The compiled library files and test cases input files included in this repository have been created with the method of using a dynamic library in mind. 
The parameters introduced at the start of all the input files can be used to modify the geometry of the specimens, as well as the material properties. No post-processing scripts have been included.
The subroutine and test cases have been run using ANSYS Mechanical APDL 2022 R2.

**INPUT**

The properties of the user-defined element, which will establish its behaviour, are introduced through the real constants assigned to each element in an ANSYS analysis. The following input must be provided for each real constant:
* 1: Number of integration points.
* 2: Flag to select the integration method: Newton-Cotes (1) or Gauss (2).
* 3 and 4: Mode I and II interface strengths (MPa).
* 5: Penalty stiffness (N/mm3).
* 6: B-K coefficient.
* 7 and 8: Mode I and II fracture toughnesses (N/mm).
* 10, 11 and 12: Paris’ law parameter C in mode I, mode II, and mixed mode (mm/cycle) (Blanco et al., Mixed-mode delamination growth in carbon-fibre composite laminates under cyclic loading. International Journal of Solids and Structures, 41:4219–4235, 7 2004. doi:10.1016/j.ijsolstr.2004.02.040).
* 13, 14 and 15: Paris’ law parameter m in mode I, mode II, and mixed mode (Blanco et al., Mixed-mode delamination growth in carbon-fibre composite laminates under cyclic loading. International Journal of Solids and Structures, 41:4219–4235, 7 2004. doi:10.1016/j.ijsolstr.2004.02.040)
* 16: Fatigue crack growth increment target (mm).
* 17 and 18: Mode I and II initiation S-N curve slopes.
* 19: Flag to indicate whether this element should undergo fatigue initiation (1) or not (0).
* 20: Length of the crack initiation region (mm)

An example of how to introduce the input can be found, for example, in the Fatigue_DNS.inp:

> *NREAL = 20 ! Number of real constants*
>
> *NSAVEVARS = (NIntP**2)*20 ! Number of saved variables*
>
> *NRSLTVAR = 0 ! Number of variables saved in results files (not including stress and total strain data)*
>
> *REAL, 2*
>
> *R, 2, NIntP, IntMethodFlag, tauI, tauII, KStiffness, BKeta*
>
> *RMORE, GIc, GIIc, visc, ParisCI, ParisCII, ParisCm*
>
> *RMORE, ParismI, ParismII, Parismm, IncTargeta, sI, sII*
>
> RMORE, 1, InitiationLength 


 
**OUTPUT**

ANSYS allows for some of the results within the subroutine to be saved as so-called "basic element results" in the 14 components of the stress and strain at each corner point of the element. This can be used to visualise the behaviour of the cohesive element by plotting these components. In the implemented user-defined cohesive element, the following results were stored in all corner points of all cohesive elements for each of the stress and strain components:
* SX. Mode I J-integral element contribution (N/mm).
* SY. Mode II J-integral element contribution (N/mm).
* SZ. Number of elements in the fracture process zone of the propagation crack region to which this element belongs in the case it does.
* SXY. Derivative of the opening separation with respect to the crack growth.
* SYZ. Derivative of the in-plane sliding separation with respect to the crack growth.
* SXZ. Traction component in the first integration point in the opening direction (MPa).
* SEQV. Traction component in the first integration point in the in-plane sliding direction (MPa).
* EPTOX. Number of cycles.
* EPTOY. Crack growth defined as total length of fully-damaged elements in the propagation crack region to which this element belongs in the case it does (mm).
* EPTOZ. Mode I energy release rate GI associated to the propagation crack region to which this element belongs in the case it does (N/mm)
* EPTOXY. Mode II energy release rate associated to the propagation crack region to which this element belongs in the case it does GI (N/mm)
* EPTOYZ. Stiffness damage variable Dk in the first integration point.
* EPTOXZ. Stiffness damage variable Dk in the last integration point.

DOCUMENTATION 
---------------
The associated journal paper (Urcelay Oca et al. An ANSYS user cohesive element for the modelling of fatigue initiation and propagation of delaminations in composite structures. Under review (2024) 2024;) contains a more in-depth description of the subroutine and the test cases.

Please send your comments or questions to Inigo Urcelay Oca: iuoc@mp.aau.dk / inigo.oca@siemensgamesa.com 


# ANSYS Fatigue Cohesive Element
This code repository contains a user element subroutine for an ANSYS Mechanical APDL cohesive element to model 2D delamination onset/initiation and propagation under quasi-static and fatigue loading. The user element calculates also the J-integral in the cohesive zone and can be used in models with multiple delaminations. Tests cases showcasing these capabilities are also included in this repository.

LICENSE: 
---------------
Copyright (c) 2023, Siemens Gamesa RE
All rights reserved.

The author would like to kindly emphasize that the following rules are respected when using the code. 
* Use of the code in source and binary forms, with or without modification are permitted provided that due reference is given to the author, the code and publications describing the implemented models.
* The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

DISCLAIMER:
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Please send your comments or questions to Iñigo Urcelay Oca: iuoc@mp.aau.dk / inigo.oca@siemensgamesa.com 

DESCRIPTION
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

DOCUMENTATION 
---------------
The subroutine and input files include comments describing the functionalities, input and output variables. The comments in the subroutine refer in some occasions to the journal papers where the implemented models are described.
<The associated journal paper ... contains a more in-depth description of the subroutine and the test cases. (Not published yet)>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! LICENSE: 
! Copyright (c) 2023, Siemens Gamesa RE
! All rights reserved.
!
! The author would like to kindly emphasize that the following rules are respected when using the code. 
!   * Use of the code in source and binary forms, with or without modification are permitted provided that due reference is given to the author, 
!   the code and publications describing the implemented models.
!   * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
!
! DISCLAIMER:
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
! Please send your comments or questions to Inigo Urcelay Oca: iuoc@mp.aau.dk / inigo.oca@siemensgamesa.com 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   This is an 8-node cohesive user-element subroutine for ANSYS to model delamination onset/initiation and propagation under quasi-static and fatigue loading. 
!   The model calculates also the J-integral in the cohesive zone, which is required for the fatigue propagation model.
!   It is based on the work described in the following papers:
!     - [1] Turon, A., Camanho, P. P., Costa, J., ; Dávila, C. G. (2006). A damage model for the simulation of delamination in advanced composites under variable-mode loading. 
!     Mechanics of Materials, 38(11), 1072–1089. https://doi.org/10.1016/j.mechmat.2005.10.003
!     - [2] Bak, B. L. V., Turon, A., Lindgaard, E., ; Lund, E. (2016). A simulation method for high-cycle fatigue-driven delamination using a cohesive zone model. 
!     International Journal for Numerical Methods in Engineering, 106(3), 163–191. https://doi.org/10.1002/nme.5117 
!     - [3] May, M., Hallett, S. R. (2010). A combined model for initiation and propagation of damage under fatigue loading for cohesive interface elements. 
!     Composites Part A: Applied Science and Manufacturing, 41(12), 1787–1796. https://doi.org/10.1016/j.compositesa.2010.08.015 
!     - [4] Camanho, P. P., Dávila, C. G., de Moura, M. F. (2003). Numerical Simulation of Mixed-mode Progressive Delamination in Composite Materials.
!     Journal of Composite Materials, 37(16), 1415–1438. https://doi.org/10.1177/002199803034505 
!     - [5] Allegri, G. (2020). A unified formulation for fatigue crack onset and growth via cohesive zone modelling. 
!     Journal of the Mechanics and Physics of Solids, 138. https://doi.org/10.1016/j.jmps.2020.103900 
!     - [6] Jiang, W. G., Hallett, S. R., Green, B. G., Wisnom, M. R. (2007). A concise interface constitutive law for analysis of delamination and splitting in composite materials and its application to scaled notched tensile specimens. 
!     International Journal for Numerical Methods in Engineering, 69(9), 1982–1995. https://doi.org/10.1002/nme.1842     
!     - [7] Tao, C., Mukhopadhyay, S., Zhang, B., Kawashita, L. F., Qiu, J., Hallett, S. R. (2018). An improved delamination fatigue cohesive interface model for complex three-dimensional multi-interface cases. 
!     Composites Part A: Applied Science and Manufacturing, 107, 633–646. https://doi.org/10.1016/j.compositesa.2018.02.008       
!   References to these papers are included throughout the subrutine. The notation differs in some occasions.
!   The static cohesive model is from [1], while the fatigue propagation model is from [2]. The initiation model is inspired by [3].
! 
!   Notes:
!     - This subroutine follows a cycle jump strategy, where the first substep is the static ramp-up up to the maximum applied force/displacement, followed by a fatigue substep, 
!     modelled on ANSYS as a substep with a constant maximum applied force/displacement
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODULE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This section creates a module where to store the data from all elements
! Some of the module variables are assigned the dimensions 1e5 and 10. For this reason, the element number for the cohesive elements must be below 1e5, and no more than 9 cracks can be modelled (since 1 is used
! for all elements not assigned to a certain crack)
!
      MODULE Module_CohElem
            IMPLICIT NONE

            INTEGER ::
     &  Substep,                                     ! Variable containing the substep number. It is used to detect when the analysis moves for the first element to a new substep
     &  Iteration,                                   ! Variable containing the substep number. It is used to detect when the analysis moves for the first element to a new iteration
     &  NumElemFPZ(10),                              ! Number of cohesive elements within the fracture process zone
     &  NumDamElem(10),                              ! Number of fully damaged cohesive elements
     &  Neighbours(1e5, 4),                          ! Variable containing the node number to the left (i,1) and to the right (i,2) of element i, and the damage status (as Dk*1E6) in the first and last integration points (i,3) and (i,4) 
     &  CrackTip(10, 2),                             ! Variable containing the element situated in the crack tip (i,1) and the direction of propagation (i,2) for each crack i
     &  CrackNumbers(1e5),                           ! Variable containing the crack number associated to each element
     &  PropCracks(10),                              ! Variable indicating whether a crack should go through fatigue propagation damage (1) or fatigue initiation damage (0)
     &  InitiatedCrack(1e5)                          ! Variable indicating whether this element belongs to an initiated crack
            DOUBLE PRECISION ::
     &  JIs(1e5, 10),                                ! Mode I J-integral element contributions
     &  JIIs(1e5, 10),                               ! Mode II J-integral element contributions
     &  cParisICrack(10),                            ! Parameter c in Paris law for mode I for this crack. Defined in case there are multiple cracks growing in various interfaces with different interface properties
     &  cParisIICrack(10),                           ! Parameter c in Paris law for mode II for this crack
     &  cParismCrack(10),                            ! Parameter c in Paris law for fitting of mode interaction for this crack
     &  mParisICrack(10),                            ! Parameter m in Paris law for mode I for this crack
     &  mParisIICrack(10),                           ! Parameter m in Paris law for mode II for this crack
     &  mParismCrack(10),                            ! Parameter m in Paris law for fitting of mode interaction for this crack
     &  GIcCrack(10),                                ! Mode I critical energy release rate / fracture toughness for this crack
     &  GIIcCrack(10),                               ! Mode II critical energy release rate / fracture toughness for this crack
     &  BKetaCrack(10),                              ! Benzeggagh–Kenane criterion coefficient for this crack
     &  JI(10),                                      ! Mode I J-integral
     &  JII(10),                                     ! Mode II J-integral
     &  CrackGrowth(10),                             ! Crack growth (based on fully damaged elements)
     &  DeltaNPropCrack(10),                         ! Increment of cycles based on the crack propagation calculated from the crack growth rate and the crack growth target per substep for each crack
     &  DeltaNInitElem(1e5),                         ! Increment of cycles based on the crack initiation for each elId
     &  SeverityInitElem(1e5),                       ! Variable containing the severity, ratio of current traction over static strength, of each element
     &  SlopeInitElem(1e5),                          ! Variable containing the initiation slope of each element based on their displacement-based mode mixity
     &  DInitElem(1e5),                              ! Variable containing the initiation damage variable of each element 
     &  DeltaNInitCrack,                             ! Increment of cycles due to initiation 
     &  DeltaN,                                      ! Increment of cycles used in the calculations based on the previous crack growth rate 
     &  ElemLengths(1e5)                             ! Variable containing the element length for every element
     
      END MODULE Module_CohElem
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START SUBROUTINE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   The user-element subroutine UserElem and its arguments are defined as done in the example of a userelem subroutine from ANSYS documentation
      subroutine UserElem (
     &                     elId, matId, keyMtx, lumpm, nDim, nNodes,
     &                     Nodes, nIntPnts, nUsrDof, kEStress, 
     &                     keyAnsMat, keySym, nKeyOpt, KeyOpt,
     &                     temper, temperB, tRef, kTherm, 
     &                     nPress, Press, kPress, nReal, RealConst, 
     &                     nSaveVars, saveVars, xRef, xCur, 
     &                     TotValDofs, IncValDofs, ItrValDofs,
     &                     VelValDofs, AccValDofs,
     &                     kfstps, nlgeom, nrkey, outkey, elPrint, iott,
     &                     keyHisUpd, ldstep, isubst, ieqitr, timval, 
     &                     keyEleErr, keyEleCnv,
     &                     eStiff, eMass, eDamp, eSStiff,
     &                     fExt, fInt, elVol, elMass, elCG, 
     &                     nRsltBsc, RsltBsc, nRsltVar, RsltVar, 
     &                     nElEng, elEnergy) 
!
!   Indicate the use of the module where the information on all cohesive elements is stored
      USE Module_CohElem  
!
!   Include the following include decks:
#include "impcom.inc" 
#include "locknm.inc"    
! 
!   The remaining lines in this section are contained in the example provided by ANSYS:
!
!     PROGRAMMER SHOULD NOT CHANGE ANY PURE INPUT ARGUMENTS (marked by ....,in)!
!
!     elId      (int,sc,in)        element number
!     matId     (int,sc,in)        material number of this element
!     keyMtx    (int,ar(10),in)    matrix and load vector form requests
!                                     0 = not requested, 1 = requested
!                                     see below for more details
!     lumpm     (int,sc,in)        mass matrix format
!                                    = 0 no lumped mass matrix
!                                    = 1 lumped mass matrix
!     nDim      (int,sc,in)        number of dimensions of the problem
!                                       (defined on USRELEM command as NDIM)
!                                    = 2 2D
!                                    = 3 3D
!     nNodes    (int,sc,in)        number of nodes of the element
!                                       (defined on USRELEM command as NNODES)
!     Nodes     (int,ar(nNodes),in)node list of this element 
!     nIntPnts  (int,sc,in)        maximum number of integration points
!                                       (defined on USRELEM command as NINTPNTS)
!     nUsrDof   (int,sc,in)        number of DOFs of this element (matrix and 
!                                     load vector size)
!     kEStress  (int,sc,in)        kEStress 
!                                       (defined on USRELEM command as KESTRESS)
!     keyAnsMat (int,sc,in)        key to indicate if ANSYS material
!                                     routine is going to be called
!                                     (defined on USRELEM command as KEYANSMAT)
!                                     = 0, No
!                                     = 1, Yes
!     keySym    (int,sc,in)        key to indicate if element matrices
!                                     is symmetric
!                                       (defined on USRELEM command as KEYSYM)
!                                     = 0, symmetric
!                                     = 1, unsymmetric
!     nKeyOpt   (int,sc,in)        number of element key options able to be
!                                     used in this routine
!     KeyOpt    (int,ar(nKeyOpt),in) values of element key option defined
!                                     by et or keyopt command for the
!                                     user elements, only the first
!                                     nKeyOpt values are passed in and can
!                                     be used to branch the routine for
!                                     different formulations
!     temper    (dp,ar(nNodes),in) nodal temperatures at current time
!     temperB   (dp,ar(nNodes),in) nodal temperatures at the beginning of this
!                                     incremental step (substep)
!     tRef      (dp,sc,in)         reference temperature
!     kTherm    (int,sc,inout)     input:  flag for thermal loading 
!                                      = 1, Temperatures at nodes are different 
!                                      from the reference temperature, 
!                                      thermal loading might be needed.
!                                      = 0, Temperatures at nodes are the same
!                                      as the reference temperature, 
!                                      thermal loading is not needed.
!                                  output:  flag for thermal strains
!     nPress    (int,sc,in)        number of pressure values for this element
!     Press     (dp,ar(nPress),in) applied elemental face load (pressure)
!     kPress    (int,sc,in)        flag for pressure loading 
!                                      = 1, pressure load is applied and 
!                                      equivalent nodal forces should be 
!                                      calculated
!                                      = 0, no pressure loading
!     nReal     (int,sc,in)        number of real constants
!                                       (defined on USRELEM command as NREAL)
!     RealConst (dp,ar(nReal),in)  user defined real constants 
!     nSaveVars (int,sc,in)        number of saved variables
!                                      (defined on USRELEM command as NSAVEVARS)
!     saveVars  (dp,ar(nSaveVars),inout) user saved variables
!     xRef      (dp,ar(nDim,nNodes),in)
!                                  nodal coordinates in initial configuration
!     xCur      (dp,ar(nDim,nNodes),in)
!                                  nodal coordinates in current configuration
!     TotValDofs (dp,ar(nUsrDof),in) total values of DOFs (displacements) 
!                                     from time = 0
!     IncValDofs (dp,ar(nUsrDof),in) incremental values of DOFs (displacements) 
!                                     for the current step
!     ItrValDofs (dp,ar(nUsrDof),in) iterative values of DOFs (displacements)
!                                     for the current iteration
!                                     (normally needed for debug only)
!     VelValDofs (dp,ar(nUsrDof),in) first time derivatives of DOFs 
!                                             (velocities) (normally not needed)
!     AccValDofs (dp,ar(nUsrDof),in) second time derivatives of DOFs 
!                                          (accelerations) (normally not needed)
!     kfstps    (int,sc,in)        key for the first iteration of first 
!                                     substep of the first load step
!                                     = 1 yes
!                                     = 0 no
!     nlgeom    (int,sc,in)        large deformation key [from nlgeom command]
!                                     = 0 NLGEOM,OFF
!                                     = 1 NLGEOM, ON
!     nrkey     (int,sc,in)        key to indicate a newton-raphson
!                                     (incremental) procedure
!                                     = 0 No
!                                     = 1 Yes
!     outkey    (int,sc,in)        key to indicate if any element output is
!                                     to be placed on the print file or the 
!                                     result file
!                                     = 0 No
!                                     = 1 Yes
!     elPrint   (int,sc,in)        key to indicate if any element output is 
!                                     to be placed on the print file
!                                     = 0 No
!                                     = 1 Yes
!     iott      (int,sc,in)        print output file unit number
!     keyHisUpd (int,sc,in)        key to indicate if history-dependent
!                                    variables need to be updated, like
!                                    equivalent plastic strain, back stress
!                                    etc. since the iteration is already
!                                    converged
!                                     = 0 not converged, don't need to update
!                                         history dependent variables
!                                     = 1 yes, converged, need to update
!                                         history dependent variables
c
!     --- The following 7 variable group can usually be ignored.
!     --- The variables are used for debug, timing, and convergence control.
!     ldstep    (int,sc,in)        current load step number
!     isubst    (int,sc,in)        current substep number
!     ieqitr    (int,sc,in)        current equilibium iteration  number
!     timval    (int,sc,in)        current time value
!     keyEleErr (int,sc,inout)     key to indicate if there is any element 
!                                     formulation error, like negative Jacobian.
!                                     The error could be caused by too
!                                     large incremental step, illegal model.
!                                     = 0 no error (preset value before calling)
!                                     = 1 some error happens. ANSYS will
!                                     decide to stop the analysis or cutback
!                                     the substep (bi-section) based on other
!                                     user input and information at higher
!                                     level.
!     keyEleCnv (int,sc,inout)     key to flag if this element satisfies
!                                     the user defined element convergence
!                                     criterion. 
!                                     = 1, yes, the criterion is satisfied
!                                       or don't have any criterion at all
!                                       it is preset value before calling
!                                     = 0, no, the element doesn't satisfy
!                                       element convergence criterion. If
!                                       this is the case, the iteration will
!                                       not converge even when both force
!                                       and displacement converge 
!       ---- end of 7 variable group -----
c
!                                                                  requested if
!     eStiff(dp,ar(nUsrDof,nUsrDof),inout) stiffness matrix         keyMtx(1)=1
!     eMass (dp,ar(nUsrDof,nUsrDof),inout) mass matrix              keyMtx(2)=1
!     eDamp (dp,ar(nUsrDof,nUsrDof),inout) damping matrix           keyMtx(3)=1
!     eSStiff(dp,ar(nUsrDof,nUsrDof),inout)stress stiffness matrix  keyMtx(4)=1
!     fExt      (dp,ar(nUsrDof),out)       applied f vector         keyMtx(5)=1
!     fInt      (dp,ar(nUsrDof),out)       internal force vector    keyMtx(6)=1

!     elVol     (dp,sc,out)        element volume
!     elMass    (dp,sc,out)        element mass
!     elCG      (dp,ar(3),out)     element centroid coordinates in current
!                                     configuration
!     nRsltBsc  (dp,sc,in)         number of basic elemental results saved in
!                                   result files 
!     RsltBsc   (dp,ar(nRsltBsc),out) basic elemental results 
!                                       (see EXPLANATION below)
!     nRsltVar  (int,sc,in)        number of elemental results saved in 
!                                     result file as non-summable miscellaneous
!                                     variables 
!                                       (defined on USRELEM command as NRSLTVAR)
!     RsltVar   (dp,ar(nRsltVar),out) variables to saved in result files as
!                                      non-summable miscellaneous variables 
!                                      requested when outkey = 1
c
!     nElEng    (int,sc,in)        number of energies (fixed at 3)
c
!     elEnergy  (dp,ar(nElEng),out) elemental energy
!                                     elEnergy(1), element strain energy
!                                     elEnergy(2), element plastic strain energy
!                                     elEnergy(3), element creep strain energy
c
!     EXPLANATION OF RsltBsc
!     
!       Basic element results are saved and total number of result 
!     quantities is nRsltBsc, where:
!            nRsltBsc = (7+7)* number of corner nodes at one element.
!       To process the quantites by post processing properly, the results 
!     must be in the following order:
!       1.) Stresses: Sx Sy Sz Sxy Syz Sxz Seqv at all corner points,
!     followed by:
!       2.) Strains : Ex Ey Ez Exy Eyz Exz Eeqv at all corner points
!     where Seqv and Eeqv = equivalent stress and strain respectively
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLE DECLARATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
!
!   Next, the USERELEM arguments are declared. Compared to ANSYS documentation: nRsltBsc is declared as integer and timval as double precision
        INTEGER ::
     &  elId,                                      ! Element number (in)
     &  matId,                                     ! Material number of this element (in)
     &  keyMtx(10),                                ! Matrix and load vector form requests (in)
     &  lumpm,                                     ! Mass matrix formatmass matrix format (in)
     &  nDim,                                      ! Number of dimensions of the problem (in)
     &  nNodes,                                    ! Number of nodes of the element (in)
     &  Nodes(nNodes),                             ! Node list of this element (in)
     &  nIntPnts,                                  ! Maximum number of integration points (in)
     &  nUsrDof,                                   ! Number of DOFs of this element (matrix and load vector size) (in)
     &  kEStress,                                  ! kEStress (in)
     &  keyAnsMat,                                 ! key to indicate if ANSYS material routine is going to be called (in)
     &  keySym,                                    ! Key to indicate if element matrices is symmetric (in)
     &  nKeyOpt,                                   ! Number of element key options able to be used in this routine (in)
     &  KeyOpt(nKeyOpt),                           ! Values of element key option defined by et or keyopt command for the user elements, only the first nKeyOpt values are passed in and can be used to branch the routine for different formulations (in)
     &  kTherm,                                    ! Flag for thermal loading (in) or strains (out)
     &  nPress,                                    ! Number of pressure values for this element (in)
     &  kPress,                                    ! Flag for pressure loading (in)
     &  nReal,                                     ! Number of real constants (in)
     &  nSaveVars,                                 ! Number of saved variables (in)
     &  kfstps,                                    ! Key for the first iteration of first substep of the first load step (in)
     &  nlgeom,                                    ! Large deformation key [from nlgeom command] (in)
     &  nrkey,                                     ! Key to indicate a Newton-Raphson (in)
     &  outkey,                                    ! Key to indicate if any element output is to be placed on the print file or the result file (in)
     &  elPrint,                                   ! Key to indicate if any element output is to be placed on the print file (in)
     &  iott,                                      ! Print output file unit number (in)
     &  keyHisUpd,                                 ! Key to indicate if history-dependent variables need to be updated, like equivalent plastic strain, back stress  etc. since the iteration is already converged (in)
     &  ldstep,                                    ! Current load step number (in)
     &  isubst,                                    ! Current substep number (in)
     &  ieqitr,                                    ! Current equilibrium iteration number (in)
     &  keyEleErr,                                 ! Key to indicate a Newton-Raphson (inout)
     &  keyEleCnv,                                 ! Key to to flag if this element satisfies the user defined element convergence criterion (inout)
     &  nRsltBsc,                                  ! Number of basic elemental results saved in result files (in)
     &  nRsltVar,                                  ! Number of elemental results saved in result file as non-summable miscellaneous variables (in)
     &  nElEng                                     ! Number of energies (fixed at 3) (in)
        DOUBLE PRECISION ::
     &  temper(nNodes),                            ! Nodal temperatures at current time (in)
     &  temperB(nNodes),                           ! Nodal temperatures at the beginning of this incremental step (substep) (in)
     &  tRef,                                      ! Reference temperature (in)
     &  Press(nPress),                             ! Applied elemental face load (pressure) (in)
     &  RealConst(nReal),                          ! User defined real constants (in)
     &  saveVars(nSaveVars),                       ! User saved variables (inout)
     &  xRef(nDim,nNodes),                         ! Nodal coordinates in initial configuration (in)
     &  xCur(nDim,nNodes),                         ! Nodal coordinates in current configuration (in)
     &  TotValDofs(nUsrDof),                       ! Total values of DOFs (displacements) from time 0 (in)
     &  IncValDofs(nUsrDof),                       ! Incremental values of DOFs (displacements) for the current step (in)
     &  ItrValDofs(nUsrDof),                       ! Iterative values of DOFs (displacements) for the current iteration (normally needed for debug only) (in)
     &  VelValDofs(nUsrDof),                       ! First time derivatives of DOFs (velocities) (normally not needed)(in)
     &  AccValDofs(nUsrDof),                       ! Second time derivatives of DOFs (accelerations) (normally not needed)(in)
     &  timval,                                    ! Current time value (in)
     &  eStiff(nUsrDof,nUsrDof),                   ! Stiffness matrix if keyMtx(1)=1 (inout)
     &  eMass(nUsrDof,nUsrDof),                    ! Mass matrix if keyMtx(2)=1 (inout)
     &  eDamp(nUsrDof,nUsrDof),                    ! Damping matrix if keyMtx(3)=1 (inout)
     &  eSStiff(nUsrDof,nUsrDof),                  ! Stress stiffness matrix if keyMtx(4)=1 (inout)
     &  fExt(nUsrDof),                             ! Applied f vector if keyMtx(5)=1 (out)
     &  fInt(nUsrDof),                             ! Internal force vector if keyMtx(6)=1 (out)
     &  elVol,                                     ! Element volume (out)
     &  elMass,                                    ! Element mass (out)
     &  elCG(3),                                   ! Element centroid coordinates in current configuration (out)
     &  RsltBsc(nRsltBsc),                         ! Basic elemental results (out) They must be in the following order: Stresses Sx Sy Sz Sxy Syz Sxz Seqv at all corner points, followed by: Strains Ex Ey Ez Exy Eyz Exz Eeqv at all corner points, where Seqv and Eeqv = equivalent stress and strain respectively
     &  RsltVar(nRsltVar),                         ! Variables to saved in result files as non-summable miscellaneous variables requested when outkey = 1 (out)
     &  elEnergy(nElEng)                           ! Elemental energy: elEnergy(1), element strain energy; elEnergy(2), element plastic strain energy; elEnergy(3), element creep strain energy (out)
!
!    Next, the variables used for debugging are declared. As well as two ANSYS subroutines regarding multithreading
        EXTERNAL ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors 
     &  pplock, ppunlock                           ! ANSYS own subroutines regarding multithreading 
      INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout                                     ! Output of erhandler (ANSYS own subroutine for displaying program errors) (inside)
!     Example of debugging using erhandler subroutine
!     The ANSYS subroutine is used to print some values in an error message. Its arguments are the filename of the overall subroutine, a number between 1-9999 
!     (recommended to use a sequence), the level of error (1 for note, 2 for warning, 3 for error, 4 for fatal, ...), the error message to display (including formatting),
!     the vector of data to display (both integer and double precision) and the vector of character data. 
!     For example, the lines below would print a warning message where the value of several variables is indicated  
!       errout = erhandler('USERELEM',1000,2,
!      & 'Printing initial values   
!      & keyMtx: %F, %F, %F, %F, %F, %F.
!      & Others: %F, %F, %F.',
!      & (/keyMtx(1)*ONE, keyMtx(2)*ONE, keyMtx(3)*ONE,  
!      & keyMtx(4)*ONE, keyMtx(5)*ONE, keyMtx(6)*ONE,
!      & nRsltBsc*ONE, isubst*ONE, ieqitr*ONE/), ' ')
!
!    Next, the variables related to the element integration are declared
        INTEGER ::
     &  NIntP,                                     ! Number of integration points (in)
     &  IntMethodFlag,                             ! Flag to select the integration method: Newton-Cotes (1) or Gauss (2) (in)
     &  IntP,                                      ! Integration point number (inside)
     &  iInt,                                      ! Integration counter along natural coordinate eta (inside)
     &  jInt                                       ! Integration counter along natural coordinate xi (inside)
        DOUBLE PRECISION ::
     &  ShapeFactors(nNodes),                      ! Shape function factors in local coordinates (inside)
     &  DerivShapeFactors(nDim-1,nNodes),          ! Factors for the derivatives of the shape functions (inside)
     &  dA,                                        ! Area scale factor for current integration point value (inside)
     &  Theta(nDim, nDim),                         ! Matrix of directional vectors. Orthogonal rotation tensor relating the local and fixed coordinate systems (inside)
     &  Bglobal(nDim, nUsrDof),                    ! B-matrix in global coordinates (inside)
     &  elJIntegralI,                              ! Mode I J-integral contribution from this element (inside)
     &  elJIntegralII                              ! Mode II J-integral contribution from this element (inside)
        DOUBLE PRECISION, ALLOCATABLE ::
     &  IntPCoords(:),                             ! Natural coordinates for the integration points (inside)
     &  IntPWeights(:),                            ! Weights for the integration points (inside)
     &  DisplJump(:,:),                            ! Displacement jump at integration point (inside)  
     &  ConstStiff(:, :, :),                       ! Constitutive stiffness matrix (inside)
     &  Traction(:, :)                             ! Traction matrix (inside)
!
!    Next, the variables related to damage are declared
        INTEGER ::
     &  NumCrack,                                  ! Crack number for this element (inside)
     &  InitiationFlag,                            ! Flag to indicate whether this element should experience fatigue initiation (1) or not (0)
     &  PropNode,                                  ! Propagation node numbering: 0 indicates node 1 (lower x coordinate) and 1 the node 2 (greater x coordinate) (inside)
     &  NeighbourElem,                             ! Neighbour element number. Used when identifying the elements belonging to each crack (inside)
     &  ExtraLeaps,                                ! Number of extra leaps used to extend the crack definition past the damaged elements (inside)
     &  IntPforInit(1)                             ! Integration point at which initiation is detected because its severity is the highest (inside)
      DOUBLE PRECISION ::
     &  tauI,                                      ! Mode I interface strength (in)
     &  tauII,                                     ! Mode II interface strength (in)
     &  KStiffness,                                ! Penalty stiffness. Initial stiffness of undamaged cohesive element (in)
     &  BKeta,                                     ! Benzeggagh–Kenane criterion coefficient (in)
     &  GIc,                                       ! Mode I critical energy release rate / fracture toughness (in)
     &  GIIc,                                      ! Mode II critical energy release rate / fracture toughness (in)
     &  cParisI,                                   ! Parameter c in Paris law for mode I (in)
     &  cParisII,                                  ! Parameter c in Paris law for mode II (in)
     &  cParism,                                   ! Parameter c in Paris law for fitting of mode interaction (in)
     &  mParisI,                                   ! Parameter m in Paris law for mode I (in)
     &  mParisII,                                  ! Parameter m in Paris law for mode II (in)
     &  mParism,                                   ! Parameter m in Paris law for fitting of mode interaction (in)
     &  GI,                                        ! Mode I energy release rate obtained from the J-integral calculation (inside)
     &  GII,                                       ! Mode II energy release rate obtained from the J-integral calculation (inside)     
     &  ddelta3da,                                 ! Derivative of the opening displacement with respect to the crack growth a (inside)
     &  ddelta1da,                                 ! Derivative of the mode II displacement with respect to the crack growth a (inside) 
     &  Deltaatarget,                              ! Target crack growth rate for the current increment of cycles (sub-in)
     &  sI,                                        ! Mode I initiation S-N curve slope (in)
     &  sII,                                       ! Mode II initiation S-N curve slope (in)
     &  InitiationLength,                          ! Length of teh crack to be introduced when initiation is detected (in)
     &  DeltaNtarget,                              ! Target increment of cycles (sub-in)   
     &  NumCycles,                                 ! Number of cycles (0 in the static-loading step) 
     &  elLength                                   ! Approximation of the element length (inside) 
        DOUBLE PRECISION, ALLOCATABLE ::
     &  KH(:),                                     ! Parameter used in the calculation of the constitutive tangent stiffness matrix (inside)
     &  DkNews(:),                                 ! Variable containing the stiffness damage variable in the current substep for all integration points (inside)
     &  SlopeInits(:),                             ! Variable containing the initiation slopes for all integration points (inside)
     &  SeverityInits(:),                          ! Variable containing the Severity, i. e., ratio of current traction over static strength (inside)
     &  DInits(:),                                 ! Variable containing the initiation damage variable for all integration points (inside)
     &  DeltaNis(:)                                ! Variable containing the cycles increment to initiation (inside)
        INTEGER, ALLOCATABLE ::   
     &  Loading(:)                                 ! Vector of flags indicating the loading or not of the cohesive element (inside)
!
!    Next, other variables are declared
       INTEGER ::
     &  i,                                         ! Used in loops (inside)
     &  j,                                         ! Used in loops (inside)
     &  k                                          ! Used in loops (inside)
!
!    Finally, declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% READING PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     
!     In order to get multi-threading to work
      CALL pplock (LOCKUPF)
!
!     Read the properties through the subroutine Properties
      CALL Properties(nReal, RealConst, NIntP, IntMethodFlag, 
     & tauI, tauII, KStiffness, BKeta, GIc, GIIc, 
     & cParisI, cParisII, cParism, mParisI, mParisII, mParism,
     & Deltaatarget, sI, sII, InitiationFlag, InitiationLength)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLE INITIATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     Zero the user saved variables and other module variables, the very first time this subroutine is run (kfstps == 1)
      IF (kfstps == 1) THEN
            saveVars(1:nSaveVars) = ZERO
            JI(:) = ZERO
            JII(:) = ZERO
            JIs(:,:) = ZERO
            JIIs(:,:) = ZERO 
            cParisICrack(:) = cParisI
            cParisIICrack(:) = cParisII
            cParismCrack(:) = cParism
            mParisICrack(:) = mParisI
            mParisIICrack(:) = mParisII
            mParismCrack(:) = mParism
            GIcCrack(:) = GIc
            GIIcCrack(:) = GIIc
            BKetaCrack(:) = BKeta
            Substep = 0
            Iteration = 0
            Neighbours(:,:) = 0
            CrackTip(:,:) = 0
            PropCracks(:) = 0
            InitiatedCrack(:) = 0
            ElemLengths(:) = ZERO
            DeltaN = ZERO
            SeverityInitElem(:) = ZERO
            SlopeInitElem(:) = ONE
            DInitElem(:) = ZERO
            ! All elements are initially assigned crack number 1 to avoid accessing vectors and matrices with 0
            CrackNumbers(:) = 1
            ! Since the minimum value will be taken, assign by default a very large value to the cycles increment
            DeltaNPropCrack(:) = 1.0D10
            DeltaNInitElem(:) = 1.0D10
            DeltaNInitCrack = 1.0D10
      END IF
!
!     In order to identify the neighbouring elements, save the node numbering of their nodes 1 and 2
      ! This is performed just once, in the second iteration of the very first substep. The very first iteration is used to zero this matrix
      IF (ldstep == 1 .AND. isubst == 1 .AND. ieqitr == 2) THEN
            Neighbours(elId, 1) = Nodes(1)
            Neighbours(elId, 2) = Nodes(2)
      END IF 
!
!     Next, for this element in particular, zero the output variables which will be then piecewise-integrated in a summation 
      ! The element tangent stiffness matrix
      DO i = 1, nUsrDof
            DO j = 1, nUsrDof
                eStiff(i, j) = ZERO
            END DO
      END DO
      ! The internal force vector
      DO i = 1, nUsrDof
            fInt(i) = ZERO
      END DO
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  NEW ITERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     When moving to a new iteration for the first element
      ! Account for the possibility that ANSYS converges in one iteration in substep n, which means Iteration won't change between substeps n and n+1. So also consider a new iteration if the iteration number has not 
      ! changed, but the substep number has
      IF (Iteration .NE. ieqitr .OR.
     & (Iteration == ieqitr .AND. Substep .NE. isubst) ) THEN
            !
            ! Call subroutine CrackIdentification to identify the elements belonging to each crack. All elements not assigned to a crack have NumCrack=1
            CALL CrackIdentification(CrackTip, CrackNumbers, 
     & Neighbours, JIs, JIIs, PropCracks, DeltaNInitElem, 
     & SeverityInitElem, DInitElem,  
     & isubst, ieqitr, elId)
            !
            ! Calculate the number of cycles for the next propagation according to the da/dN based on the value of the J-integral from the previous substep through the Paris' law
            ! The subroutine CyclesPropagation is called for this calculation for each of the cracks
            ! This is only done in the first iteration so that only the J-integral from the converged solution is employed. keyHisUpd is used to avoid the possible converged iteration numbered as 1 
            IF (ieqitr == 1 .AND. keyHisUpd .NE. 1) THEN
                  ! Loop through the cracks
                  DO i = 2, SIZE(CrackTip, DIM = 1)
                        ! Calculate the number of cycles only for propagation cracks
                        IF (PropCracks(i) .NE. 0) THEN
                              CALL CyclesPropagation(
     & SUM(JIs(:, i)), SUM(JIIs(:, i)), 
     & cParisICrack(i), cParisIICrack(i), cParismCrack(i), 
     & mParisICrack(i), mParisIICrack(i), mParismCrack(i),  
     & GIcCrack(i), GIIcCrack(i), BKetaCrack(i), 
     & Deltaatarget, DeltaNPropCrack(i),
     & isubst, ieqitr, elId)
                        END IF
                  END DO
            END IF
            !
            ! Calculate the number of cycles in this increment
            ! Only in the fatigue substep
            IF (ldstep > 1) THEN
                  ! Call subroutine AdjustCyclesPropagation to adjust the cycle increment according to propagation in the crack with the highest da/dN. This subroutine modifies the cycle increment DeltaN
                  ! to aid convergence, depending on the number of iterations
                  ! If the iteration has converged, use the iteration number from the previous iteration. This is done to avoid changing in the converged iteration the input of the number of iterations, because
                  ! this number decreases by 1 in the converged iteration with respect to the previous iteration. 
                  CALL AdjustCyclesPropagation(   
     & MINVAL(DeltaNPropCrack(2:10)), DeltaN, 
     & isubst, ieqitr+1*keyHisUpd, elId)
                  !
                  ! Use the initiation cycle increment if it is lower than the one calculated for propagation
                  ! This is only done after the first iteration, since DeltaNInitCrack is calculated at the start of each substep a few lines below
                  ! Use keyHisUpd to only avoid the iteration numbered as 1 which has not converged. That is, only run the followng lines if the iteration is not 1 or it is the converged one
                  IF (ieqitr .NE. 1 .OR. keyHisUpd == 1) THEN
                        IF (DeltaN > DeltaNInitCrack) THEN
                              DeltaN = DeltaNInitCrack
                        END IF
                  END IF
            END IF
            !
            ! Save the current iteration number, so that for the rest of elements of this type, these previous lines are not run
            Iteration = ieqitr
      END IF
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  NEW SUBSTEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     When moving to a new substep for the first element
      IF (Substep .NE. isubst) THEN
            !
            ! Calculate for the defined cracks the mode I and II J-integral components
            DO i = 1, SIZE(CrackTip, DIM = 1)
                  ! For the crack numbers with elements associated to them
                  IF (CrackTip(i,1) .NE. 0) THEN
                        !
                        ! Sum up the J-integral element contributions for each crack
                        JI(i) = SUM(JIs(:, i)) 
                        JII(i) = SUM(JIIs(:, i)) 
                  END IF
            END DO
            !
            ! In the fatigue substep, study whether fatigue initiation takes place through subroutine FatigueInitiation 
            ! Only call this subroutine if at least the minimum number of cycles to initiation for all elements is lower than the minimum cycle increment according to propagation            
            IF (ldstep > 1 .AND. 
     & MINVAL(DeltaNInitElem) < MINVAL(DeltaNPropCrack(2:10))) THEN
                  ! Call subroutine FatigueInitiation to determine whether initiation takes place
                  ! This subroutine first calculates the DeltaN from initiation and then compares it to the one from propagation, making DeltaN equal to the lowest of the two
                  ! If initiation takes place, this subroutine defines the elements belonging to the 'initiation crack'. 
                  CALL FatigueInitiation(InitiationLength, ElemLengths, 
     & Neighbours, InitiatedCrack, DeltaNInitElem, DeltaNInitCrack,  
     & SeverityInitElem, SlopeInitElem, DInitElem, DeltaN,  
     & isubst, ieqitr, elId)
            ! If initiation does not take place, replace the possible low DeltaNInitCrack from the previous substep by an artificially high value
            ELSEIF (MINVAL(DeltaNInitElem) > 
     & MINVAL(DeltaNPropCrack(2:10))) THEN
                  DeltaNInitCrack = 1.0D10
            END IF
            !
            ! Zero the variables updated in every substep depending on the converged damage status of the cohesive elements 
            ! The variable containing the number of elements in the fracture process zone (FPZ)
            NumElemFPZ(:) = 0
            ! And the fully damaged ones
            NumDamElem(:) = 0
            ! And the crack growth
            CrackGrowth(:) = 0
            !
            ! Save the current subset number, so that for the rest of elements of this type, these previous lines are not run
            Substep = isubst
      END IF
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLE ALLOCATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     The variables whose dimensions depend on the number of integrated points are allocated in this section
      ALLOCATE(IntPCoords(NIntP))                  ! Natural coordinates for the integration points (inside)
      ALLOCATE(IntPWeights(NIntP))                 ! Factors for the derivatives of the shape functions (inside)
      ALLOCATE(DisplJump(NIntP**2, nDim))          ! Displacement jump at the integration points (inside)
      ALLOCATE(ConstStiff(NIntP**2, nDim, nDim))   ! Constitutive stiffness matrix (inside)
      ALLOCATE(KH(NIntP**2))                       ! Parameter used in the calculation of the constitutive tangent stiffness matrix (inside)
      ALLOCATE(Loading(NIntP**2))                  ! Vector of flags indicating the loading or not of the cohesive element (inside)
      ALLOCATE(Traction(NIntP**2, nDim))           ! Traction matrix (inside)
      ALLOCATE(DkNews(NIntP**2))                   ! Variable containing the stiffness damage variable in the current substep for all integration points (inside)
      ALLOCATE(SlopeInits(NIntP**2))               ! Variable containing the initiation slopes for all integration points (inside)
      ALLOCATE(SeverityInits(NIntP**2))            ! Variable containing the Severity, i. e., ratio of current traction over static strength (inside)
      ALLOCATE(DInits(NIntP**2))                   ! Variable containing the initiation damage variable for all integration points (inside)
      ALLOCATE(DeltaNis(NIntP**2))                 ! Variable containing the cycles increment to initiation (inside)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN BODY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     Determine the crack number for this element 
      NumCrack = CrackNumbers(elId)
      ! If the element cannot go through fatigue initiation, define this crack as a propagation one
      IF (InitiationFlag == 0) THEN
            PropCracks(NumCrack) = 1
      END IF
!
!     Obtain for this delamination the mode I and mode II energy release rates from the J-integral calculation. These are employed for damage calculations under fatigue loading
!     Only do this for elements belonging to a crack. If not, assign some fixed values just to avoid causing any issues in the intermediate calculations even if the values wouldn't affect the final results
      IF (NumCrack .NE. 1) THEN
            GI = JI(NumCrack) 
            GII = JII(NumCrack) 
      ELSE
            GI = 1.0D-10
            GII = ZERO
      END IF   
!
!     Calculate the derivatives of the opening and shearing displacements with respect to the crack growth (ddelta3da & ddelta1da) through subroutine Slopes
!     These derivatives are used later in the calculation of the fatigue damage and of the J-integral
      ! Read first the node indicating the propagation direction. A value of 0 is used if the element does not belong to a crack
      IF (NumCrack .NE. 1) THEN
            PropNode = CrackTip(CrackNumbers(CrackTip(NumCrack,1)),2)
      ELSE 
            PropNode = 0
      END IF 
      CALL Slopes(xRef, nDim, nNodes, nlgeom, xCur, 
     & TotValDofs, nUsrDof, PropNode, ddelta3da, ddelta1da, 
     & isubst, ieqitr, elId)
!
!     Obtain the integration points (coordinates and weights) through subroutine IntegrationPoints      
      CALL IntegrationPoints(NIntP, IntPCoords, IntPWeights, 
     & IntMethodFlag)
!
!     Loop over the integration points natural coordinates to calculate the element tangent stiffness matrix and the internal force vector through a piecewise-integration
      IntP = 0
      ! Loop along the element natural coordinate eta
      DO iInt = 1, NIntP
         ! Loop along the element natural coordinate xi
         DO jInt = 1, NIntP
            IntP = IntP + 1
            ! Evaluate the shape functions and its derivatives for each integration point through subroutine ShapeFunc
            CALL ShapeFunc(IntPCoords(iInt), IntPCoords(jInt), 
     &  ShapeFactors, DerivShapeFactors, nNodes, nDim)
            ! Calculate the Jacobian matrix and obtain as well the area scale factor and the orthogonal rotation tensor relating the local and fixed coordinate systems for the current integration point 
            ! through subroutine Jacobian
            CALL Jacobian(nDim, nNodes, DerivShapeFactors, xRef, nlgeom,
     &  xCur, Theta, dA, NIntP, IntPWeights, iInt, jInt)
            ! Calculate the displacements in local coordinates (and the B matrix in global coordinates) through subroutine Displacements
            CALL Displacements(nDim, nUsrDof, nNodes, ShapeFactors,   
     &  Bglobal, Theta, DisplJump, NIntP, IntP, TotValDofs)
            ! Calculate the damage variable and thus the stiffness constitutive matrix and parameter KH for the calculation of the tangent matrix through subroutine Damage
            CALL Damage(nDim, NIntP, IntP, DisplJump, tauI, tauII,  
     &  KStiffness, BKeta, GIc, GIIc, nSaveVars, saveVars, ldstep,    
     &  Loading, DkNews(IntP), timval, GI, GII, PropCracks(NumCrack), 
     &  cParisI, cParisII, cParism, mParisI, mParisII, mParism, DeltaN, 
     &  kfstps, ddelta3da, ddelta1da, SlopeInits(IntP),    
     &  SeverityInits(IntP), DInits(IntP), sI, sII, DeltaNis(IntP), 
     &  InitiatedCrack(elId), ConstStiff, KH, keyHisUpd,   
     &  ieqitr, isubst, elId)
            ! Calculate the element tangent stiffness matrix through subroutine TangentStiffness
            CALL TangentStiffness(nDim, NIntP, IntP, nUsrDof, 
     & Loading, ConstStiff, DisplJump, KH, eStiff, Bglobal, dA)
            ! Calculate the internal force vector in global coordinates through subroutine InternalForce
            CALL InternalForce(nDim, Traction, NIntP, IntP,  
     & nUsrDof, fInt, ConstStiff, DisplJump, Bglobal, dA)
         END DO
      END DO
!
!     Save the damage status from the first and last integration points in module variable Neighbours for each element.
      ! Since Neighbours mis defined to conain only integers, the value of Dk is stored multiplied by 1E6, so that the specific value of Dk can also be known. 
      Neighbours(elId, 3) = INT(DkNews(1)*1E6)
      Neighbours(elId, 4) = INT(DkNews(NIntP**2)*1E6) 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONVERGED ITERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     If the iteration has converged
      IF (keyHisUpd == 1) THEN
            !
            ! Calculate the J-integral element contributions through subroutine JIntegralCalc
            CALL JIntegralCalc(nDim, nNodes, xCur, NIntP,
     & elJIntegralI, elJIntegralII, elLength, ddelta3da, 
     & ddelta1da, Traction, IntPWeights, elId, isubst) 
            ! Save the element contributions in the module variables JIs and JIIs, 
            JIs(elId, NumCrack) = elJIntegralI
            JIIs(elId, NumCrack) = elJIntegralII
            ! Save the value of the element length for each element in the module variable ElemLengths
            ElemLengths(elId) = elLength
            !
            ! Sum the number of elements in the fracture process zone FPZ, i.e., those elements where the damage variables is > 0 and < 1   
            DO i = 1, NIntP**2
                  ! The elements are identified through the stiffness damage variable DkNew, saved in saveVars(IntP)
                  ! If the Dk of any integration point is > 0 and < 1, then the element is considered as part of the FPZ
                  IF (saveVars(i) < ONE .AND. 
     & saveVars(i) > 1.0D-6) THEN
                        NumElemFPZ(NumCrack) = 
     & NumElemFPZ(NumCrack) + 1 
                        EXIT
                  END IF
            END DO
            !
            ! Modify the variables of the elements belonging to the crack 
            ! If the element is not part of a propagation crack, save the values required to calculate the cycle increment to initiation from those obtained for all integration points to one sole per element
            IF (PropCracks(NumCrack) .NE. 1) THEN
                  ! Save the number of cycles to initiation of each element as the minimum of all the integration points
                  DeltaNInitElem(elId) = MINVAL(DeltaNis)
                  ! Save the severity of each element as the maximum of all the integration points
                  SeverityInitElem(elId) = MAXVAL(SeverityInits)
                  ! Save the initiation slope of each element as the one for the integration point where the severity is the highest
                  IntPforInit = MAXLOC(SeverityInits)
                  SlopeInitElem(elId) = SlopeInits(IntPforInit(1))
                  ! And save the initiation damage variable of each element as the maximum value
                  DInitElem(elId) = MAXVAL(DInits)
            ! If the element belongs to a propagation crack, avoid any initiation occuring in this element by saving an artificially high number of cycles to initiation
            ELSE
                  DeltaNInitElem(elId) = 1.0D10
                  ! Modify as well the other values used to calculate DeltaNInit 
                  SeverityInitElem(elId) = ZERO
                  DInitElem(elId) = -1.0D10
            END IF
            !
            ! Sum to the previous number of cycles, stored in saveVars, the used cycles increment, and store the new total in saveVars
            NumCycles = saveVars((NIntP**2)*10 + 1)
            NumCycles = NumCycles + DeltaN
            saveVars((NIntP**2)*10 + 1) = NumCycles
            !
            ! Calculate the current crack growth based on the elements fully damaged
            ! The elements are understood to be fully-damaged when all integration points have DkNew=1. So the sum of their DkNews must be greater than 0.99999999d0*NIntP**2
            IF (SUM(saveVars(1:NIntP**2)) > 
     & 0.99999999d0*NIntP**2) THEN
                  NumDamElem(NumCrack) = NumDamElem(NumCrack) + 1 
                  CrackGrowth(NumCrack) = CrackGrowth(NumCrack) 
     & + elLength
                  ! If the elements are fully damaged, zero as well their J-integral contributions
                  ! It can be the case that they aren't directly zero if Displ(IntP,3) is negative and thus the opening tractions are not zero even if the elements are fully damaged
                  JIs(elId, NumCrack) = ZERO
                  JIIs(elId, NumCrack) = ZERO
            END IF
            !
            ! Detect if the current element is the one at a crack tip, i.e, the damage variable Dk is 0 in one integration point and between 1 and 0 in another one on the other side
            CALL CrackTipIdentification(NIntP, DkNews, CrackNumbers,
     & CrackTip, elId, cParisICrack, cParisI, cParisIICrack, 
     & cParisII, cParismCrack, cParism, mParisICrack, mParisI,
     & mParisIICrack, mParisII, mParismCrack, mParism, GIcCrack, 
     & GIc, GIIcCrack, GIIc, BKetaCrack, BKeta,
     & isubst, ieqitr)
      END IF
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR POSTPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     Save internal variables in the result basic variables variables meant for stresses and strains
      ! These variables can then be plotted on ANSYS by referring to the specific stress or strain component
      ! Loop through the 8 nodes, saving the same elemnet values in all of them
      DO i = 1, 8
            ! Stresses (Sx Sy Sz Sxy Syz Sxz Seqv)
            RsltBsc(1+7*(i-1)) = elJIntegralI ! SX
            RsltBsc(2+7*(i-1)) = elJIntegralII ! SY
            RsltBsc(3+7*(i-1)) = NumElemFPZ(NumCrack) ! SZ
            RsltBsc(4+7*(i-1)) = ddelta3da ! SXY           
            RsltBsc(5+7*(i-1)) = ddelta1da ! SYZ
            RsltBsc(6+7*(i-1)) = Traction(1,3) ! SXZ           
            RsltBsc(7+7*(i-1)) = Traction(1,1) ! SEQV
            ! Strains(Ex Ey Ez Exy Eyz Exz Eeqv)
            RsltBsc(1+56+7*(i-1)) = NumCycles ! EPTOX
            RsltBsc(2+56+7*(i-1)) = CrackGrowth(NumCrack) ! EPTOY
            RsltBsc(3+56+7*(i-1)) = GI ! EPTOZ
            RsltBsc(4+56+7*(i-1)) = GII ! EPTOXY           
            RsltBsc(5+56+7*(i-1)) = NumCrack ! EPTOYZ
            RsltBsc(6+56+7*(i-1)) = DkNews(1) ! EPTOXZ           
            RsltBsc(7+56+7*(i-1)) = DkNews(NIntP**2) ! EPTOEQV
      END DO
!
      CALL ppunlock (LOCKUPF)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USERELEM SUBROUTINE END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      RETURN
      END
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     Subroutine to read the properties
      SUBROUTINE Properties(nReal, RealConst, NIntP, IntMethodFlag, 
     & tauI, tauII, KStiffness, BKeta, GIc, GIIc, 
     & cParisI, cParisII, cParism, mParisI, mParisII, mParism,
     & Deltaatarget, sI, sII, InitiationFlag, InitiationLength)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      INTEGER 
     &  nReal,                                     ! Number of real constants (sub-in)
     &  NIntP,                                     ! Number of integration points (sub-out)
     &  IntMethodFlag,                             ! Flag to select the integration method: Newton-Cotes (1) or Gauss (2) (sub-out)
     &  InitiationFlag                             ! Flag to indicate whether this element should experience fatigue initiation (1) or not (0)
      DOUBLE PRECISION
     &  RealConst(nReal),                          ! User defined real constants (sub-in)
     &  tauI,                                      ! Mode I interface strength (sub-out)
     &  tauII,                                     ! Mode II interface strength (sub-out)
     &  KStiffness,                                ! Penalty stiffness. Initial stiffness of undamaged cohesive element (sub-out)
     &  BKeta,                                     ! Benzeggagh–Kenane criterion coefficient (sub-out)
     &  GIc,                                       ! Mode I critical energy release rate / fracture toughness (sub-out)
     &  GIIc,                                      ! Mode II critical energy release rate / fracture toughness (sub-out)
     &  cParisI,                                   ! Parameter c in Paris law for mode I (sub-out)
     &  cParisII,                                  ! Parameter c in Paris law for mode II (sub-out)
     &  cParism,                                   ! Parameter c in Paris law for fitting of mode interaction (sub-out)
     &  mParisI,                                   ! Parameter m in Paris law for mode I (sub-out)
     &  mParisII,                                  ! Parameter m in Paris law for mode II (sub-out)
     &  mParism,                                   ! Parameter m in Paris law for fitting of mode interaction (sub-out)
     &  Deltaatarget,                              ! Target crack growth rate for the current increment of cycles (sub-out)
     &  sI,                                        ! Mode I initiation S-N curve slope (sub-out)
     &  sII,                                       ! Mode II initiation S-N curve slope (sub-out)
     &  InitiationLength                           ! Length of the crack to be introduced when initiation is detected (sub-out)
!    Next, the variables used in debugging are declared
        EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
        INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout                                     ! Output of erhandler (ANSYS own subroutine for displaying program errors) (inside)
!
!     Read the properties, which have been introduced as real constants
      ! Related to the integration of the cohesive element
      NIntP = IDINT(RealConst(1))                                     
      IntMethodFlag = IDINT(RealConst(2))  
      ! Related to damage under static loading
      tauI = RealConst(3)                                        
      tauII = RealConst(4)                                       
      KStiffness = RealConst(5)                                  
      BKeta = RealConst(6)                                 
      GIc = RealConst(7)                                      
      GIIc = RealConst(8)
      ! Related to damage under fatigue loading
      cParisI = RealConst(10)
      cParisII = RealConst(11)
      cParism = RealConst(12)
      mParisI = RealConst(13)
      mParisII = RealConst(14)
      mParism = RealConst(15)
      Deltaatarget = RealConst(16)
      sI = RealConst(17)
      sII = RealConst(18)
      ! Indicating whether the element should experience fatigue initiation or not
      InitiationFlag = INT(RealConst(19))
      ! The length of the fatigue initiation crack
      InitiationLength = RealConst(20)
!
!     End of the subroutine
      RETURN
      END SUBROUTINE Properties
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE CrackIdentification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     Subroutine to identify the elements belonging to each crack. All elements not assigned to a crack have NumCrack=1
      SUBROUTINE CrackIdentification(CrackTip, CrackNumbers, 
     & Neighbours, JIs, JIIs, PropCracks, DeltaNInitElem, 
     & SeverityInitElem, DInitElem,  
     & isubst, ieqitr, elId)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      INTEGER     
     &  CrackTip(10, 2),                           ! Variable containing the element situated in the crack tip (i,1) and the direction of propagation (i,2) for each crack i (sub-in)
     &  CrackNumbers(1e5),                         ! Variable containing the crack number associated to each element (sub-inout)
     &  PropNode,                                  ! Propagation node numbering: 0 indicates node 1 (lower x coordinate) and 1 the node 2 (greater x coordinate) (inside)
     &  NeighbourElem,                             ! Neighbour element number. Used when identifying the elements belonging to each crack (inside)
     &  Neighbours(1e5, 4),                        ! Variable containing the nodes to left (i,1) and right (i,2) of element i, and the damage status in the last integration point (i,3) (sub-in)
     &  PropCracks(10),                            ! Variable indicating whether a crack should go through fatigue propagation damage (1) or fatigue initiation damage (0) (sub-out)     
     &  ExtraLeaps,                                ! Number of extra leaps in the loop to identify neighbouring elements. This is used to extend the crack definition past the damaged elements (inside)
     &  i,                                         ! Used in loops (inside)
     &  j,                                         ! Used in loops (inside)
     &  k                                          ! Used in loops (inside)
      DOUBLE PRECISION
     &  JIs(1e5, 10),                              ! Mode I J-integral element contributions (sub-inout)
     &  JIIs(1e5, 10),                             ! Mode II J-integral element contributions (sub-inout)
     &  DeltaNInitElem(1e5),                       ! Increment of cycles based on the crack initiation for each elId (sub-inout)
     &  SeverityInitElem(1e5),                     ! Variable containing the severity, ratio of current traction over static strength, of each element (sub-inout)
     &  DInitElem(1e5)                             ! Variable containing the initiation damage variable of each element (sub-inout)
!    Next, the variables used in debugging are declared
      EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
      INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout,                                    ! Output of erhandler (ANSYS own subroutine for displaying program errors)                                 
     &  isubst,                                    ! Current substep number (sub-in)
     &  ieqitr,                                    ! Current equilibrium iteration number (in)
     &  elId                                       ! Element number (sub-in)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
      ! Loop through the elements at the crack tips and identify the elements in each crack contributing to the J-integral
      ! The crack tip element is defined as that where the integration points on one side are damaged, while not on the opposite side
      DO i = 2, SIZE(CrackTip, DIM = 1)
            ! If there is not an element number saved for crack number i as being the crack tip, and this element has not yet been asigned a crack number, the crack is being identified for the first time
            ! Nevertheless, this part is still run in the case that it is not a propagation crack (to avoid not detecting propagation cracks emerging from initiation cracks if the new crack tip belongs to an assigned crack)
            IF (CrackTip(i,1) .NE. 0 .AND. 
     & (CrackNumbers(CrackTip(i,1)) .NE. i .OR. 
     & PropCracks(i) .NE. 1)) THEN 
                  ! Assign the crack number to the element at the crack tip
                  CrackNumbers(CrackTip(i,1)) = i
                  !
                  ! Identify the elements of the crack tip in the damaged region
                  ! Read the node position indicating the propagation direction and employ the opposite direction
                  IF (CrackTip(i,2) == 1) THEN
                        PropNode = 1
                  ELSE
                        PropNode = 0
                  END IF
                  ! Save the element J-integral contribution of the element at the crack tip for this crack number (i)
                  ! Modify the sign if necessary because since the elements hadn't yet been assigned to a crack, the propagation direction used in the calculation of the contributions may not have been the right one
                  ! This is done through (TWO*PropNode-1), which is 1 for PropNode=1 and -1 if PropNode=0
                  JIs(CrackTip(i,1), i) = (-TWO*PropNode+ONE)*
     & JIs(CrackTip(i,1), CrackNumbers(CrackTip(i,1)))
                  JIIs(CrackTip(i,1), i) = (-TWO*PropNode+ONE)*
     & JIIs(CrackTip(i,1), CrackNumbers(CrackTip(i,1)))
                  ! Set the crack tip element as the first element from which to identify all the other neighbours
                  NeighbourElem = CrackTip(i,1)
                  ! Loop through the elements to find the neighbouring damaged elements in the direction opposite to propagation
                  LoopDamNeighs: DO k = 1, 1e5
                        LoopOneDamNeigh: DO j = 1, 1e5
                              IF (Neighbours(j, PropNode+1) == 
     & Neighbours(NeighbourElem, ABS(PropNode-1)+1)) THEN
                                    ! Check now that the element is partly damaged and more than the previous (in case damage does not initiate at an edge)
                                    ! Depending on the propagation direction, use the damage data from the integration point at one side or the other
                                    IF ((PropNode == 0 .AND.
     & Neighbours(j,3) .GE. 1 .AND. Neighbours(j,3) < 1E6 .AND.
     & Neighbours(j,3) .GE. Neighbours(NeighbourElem,3)) 
     & .OR. 
     & (PropNode == 1 .AND.
     & Neighbours(j,4) .GE. 1 .AND. Neighbours(j,4) < 1E6 .AND.
     & Neighbours(j,4) .GE. Neighbours(NeighbourElem,4))
     & ) THEN
                                          ! Save the J-integral contributions of the element (j) for this crack number (i)
                                          JIs(j, i) = 
     & (-TWO*PropNode+ONE)*JIs(j, CrackNumbers(j))
                                          JIIs(j, i) = 
     & (-TWO*PropNode+ONE)*JIIs(j, CrackNumbers(j))
                                          ! Assign the crack number (i) to the element (j)
                                          CrackNumbers(j) = i
                                          ! Set the new element as the new NeighbourElem to be used in the next iteration 
                                          NeighbourElem = j 
                                          EXIT LoopOneDamNeigh
                                    ! The next neighbour element is less damaged. This is the case for cracks not originating at an edge 
                                    ! This could be because two cracks are generated from a stress concentration far from the edge. Stop then before identifying both crackas as one sole
                                    ! Or because it is a propagation crack originating from an initiation crack
                                    ELSE
                                          ! If the loop is exited because the next neighbour is fully damaged, then this is a fatigue damage propagation crack
                                          IF (
     & Neighbours(j,3) > 1E6-1.0D-10 .AND.
     & Neighbours(j,4) > 1E6-1.0D-10) THEN
                                             PropCracks(i) = 1
                                          END IF
                                          ! Exit the loop 
                                          EXIT LoopDamNeighs
                                    END IF 
                              END IF
                        END DO LoopOneDamNeigh
                        ! Exit the loop if there are no more neighbouring elements, i. e., the iteration through all elements hasn't found a match
                        IF (j == 1e5+1) THEN
                              EXIT LoopDamNeighs
                        END IF
                  END DO LoopDamNeighs
                  !
                  ! Find the subsequent neighbouring elements in the propagation direction to correctly calculate for the J-integral (their contribution will generally be very small)
                  NeighbourElem = CrackTip(i,1)
                  ! Zero the number of leaps, which is used to extend the crack a certain number of elements past the damaged ones
                  ExtraLeaps = 0
                  ! Read the node position indicating the propagation direction 
                  IF (CrackTip(i,2) == 1) THEN 
                        PropNode = 0
                  ELSE
                        PropNode = 1
                  END IF
                  ! Loop through the elements to find the neighbouring damaged elements in the propagation direction
                  LoopPropNeighs: DO k = 1, 1e5
                        LoopOnePropNeigh: DO j = 1, 1e5
                              IF (Neighbours(j, PropNode+1) == 
     & Neighbours(NeighbourElem, ABS(PropNode-1)+1)) THEN
                                    ! Exit the loop if the next element is more damaged than the previous, such as when two cracks coincide
                                    ! Depending on the propagation direction, use the damage data from the integration point at one side or the other
                                    IF ((PropNode == 0 .AND.
     & ((Neighbours(j,3) - Neighbours(NeighbourElem,3)) .GE. 1.0D-5))
     & .OR. 
     & (PropNode == 1 .AND.
     & ((Neighbours(j,4) - Neighbours(NeighbourElem,4)) .GE. 1.0D-5))
     & ) THEN
                                          EXIT LoopPropNeighs
                                    END IF
                                    ! Exit the loop after a certain number of leaps if the next neighbour is not damaged
                                    ! Check that the number of leaps has not been exceeded
                                    IF (ExtraLeaps < 400) THEN                                          
                                          ! Save the J-integral contribution and crack number
                                          ! Since PropNode has now 
                                          JIs(j, i) = 
     & (TWO*PropNode-ONE)*JIs(j, CrackNumbers(j))
                                          JIIs(j, i) = 
     & (TWO*PropNode-ONE)*JIIs(j, CrackNumbers(j))
                                          ! Save this element number (j) as the neighbour element
                                          CrackNumbers(j) = i
                                          NeighbourElem = j 
                                          ExtraLeaps = ExtraLeaps + 1
                                          EXIT LoopOnePropNeigh
                                    ELSE
                                          EXIT LoopPropNeighs
                                    END IF
                              END IF
                        END DO LoopOnePropNeigh
                        ! Exit the loop if there are no more neighbouring elements
                        IF (j == 1e5+1) THEN
                              EXIT LoopPropNeighs
                        END IF
                  END DO LoopPropNeighs
            ! If the element has already been assigned to a crack, the crack has already been identified before. Extend its definition from the new crack tip
            ! Only for propagation cracks (to avoid not detecting propagation cracks emerging from initiation cracks if the new crack tip belongs to an assigned crack)
            ELSE IF (CrackTip(i,1) .NE. 0 .AND. 
     & CrackNumbers(CrackTip(i,1)) == i .AND. 
     & PropCracks(i) == 1) THEN
                  ! Find the subsequent neighbouring elements in the propagation direction
                  NeighbourElem = CrackTip(i,1)
                  ! Zero the number of leaps, which is used to extend the crack a certain number of elements past the damaged ones
                  ExtraLeaps = 0
                  ! Read the node position indicating the propagation direction 
                  IF (CrackTip(i,2) == 1) THEN 
                        PropNode = 0
                  ELSE
                        PropNode = 1
                  END IF
                  LoopPropNeighs2: DO k = 1, 1e5
                        LoopOnePropNeigh2: DO j = 1, 1e5
                              IF (Neighbours(j, PropNode+1) == 
     & Neighbours(NeighbourElem, ABS(PropNode-1)+1)) THEN
                                    ! Exit the loop if the next element is more damaged than the previous, such as when two cracks coincide
                                    ! Depending on the propagation direction, use the damage data from the integration point at one side or the other
                                    IF ((PropNode == 0 .AND.
     & ((Neighbours(j,3) - Neighbours(NeighbourElem,3)) .GE. 1.0D-5)) 
     & .OR. 
     & (PropNode == 1 .AND.
     & ((Neighbours(j,4) - Neighbours(NeighbourElem,4)) .GE. 1.0D-5))
     & ) THEN
                                          EXIT LoopPropNeighs2
                                    END IF
                                    ! Exit the loop after a certain number of leaps if the next neighbour is not damaged
                                    ! Check that the number of leaps has not been exceeded
                                    IF (ExtraLeaps < 400) THEN                                          
                                          ! Save the J-integral contribution and crack number
                                          JIs(j, i) = 
     & JIs(j, CrackNumbers(j))
                                          JIIs(j, i) = 
     & JIIs(j, CrackNumbers(j))
                                          ! Save this element number (j) as the neighbour element
                                          CrackNumbers(j) = i
                                          NeighbourElem = j 
                                          ExtraLeaps = ExtraLeaps + 1
                                          EXIT LoopOnePropNeigh2
                                    ELSE
                                          EXIT LoopPropNeighs2
                                    END IF
                              END IF
                        END DO LoopOnePropNeigh2
                        ! Exit the loop if there are no more neighbouring elements
                        IF (j == 1e5+1) THEN
                              EXIT LoopPropNeighs2
                        END IF
                  END DO LoopPropNeighs2
            END IF
      END DO
      !
!     Modify some fatigue initiation variables to artificial values if the elements are identified as being part of a propagation crack
      ! This is done for those propagation cracks generated from initiation cracks, such that initiation won't be activated in them
      DO i = 1, 1e5
            IF (PropCracks(CrackNumbers(i)) == 1) THEN
                  DeltaNInitElem(i) = 1.0D10
                  SeverityInitElem(i) = ZERO
                  DInitElem(i) = -1.0D10
            END IF
      END DO
!
!     End of the subroutine
      RETURN
      END SUBROUTINE CrackIdentification
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE CyclesPropagation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Subroutine to calculate the cycle increment for each propagation crack
      SUBROUTINE CyclesPropagation(GI, GII, cParisI, cParisII, 
     & cParism, mParisI, mParisII, mParism, GIc, GIIc, BKeta, 
     & Deltaatarget, DeltaNProp, isubst, ieqitr, elId)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      DOUBLE PRECISION
     &  G,                                         ! Total energy release rate obtained from the J-integral calculation (inside)
     &  GI,                                        ! Mode I energy release rate obtained from the J-integral calculation (sub-in)
     &  GII,                                       ! Mode II energy release rate obtained from the J-integral calculation (sub-in)
     &  modemixity,                                ! Mode mixity based on the J-integral calculation (inside)
     &  cParis,                                    ! Parameter c in Paris law for a certain mode mixity (inside)
     &  cParisI,                                   ! Parameter c in Paris law for mode I (sub-in)
     &  cParisII,                                  ! Parameter c in Paris law for mode II (sub-in)
     &  cParism,                                   ! Parameter c in Paris law for fitting of mode interaction (sub-in)
     &  mParis,                                    ! Parameter m in Paris law for a certain mode mixity (inside)
     &  mParisI,                                   ! Parameter m in Paris law for mode I (sub-in)
     &  mParisII,                                  ! Parameter m in Paris law for mode II (sub-in)
     &  mParism,                                   ! Parameter m in Paris law for fitting of mode interaction (sub-in)
     &  GcGlobal,                                  ! Mixed mode global critical energy release rate (inside)    
     &  GIc,                                       ! Mode I critical energy release rate / fracture toughness (sub-in)
     &  GIIc,                                      ! Mode II critical energy release rate / fracture toughness (sub-in)
     &  BKeta,                                     ! Benzeggagh–Kenane criterion coefficient (sub-in)
     &  dadN,                                      ! Crack growth rate: derivative of the crack growth wrt the number of cycles (inside)
     &  Deltaatarget,                              ! Target crack growth rate for the current increment of cycles (sub-in)
     &  DeltaNProp                                 ! Cycles increment according to the current crack growth rate and crack growth target per substep (sub-out)
!    Next, the variables used in debugging are declared. Including some variables used to identify the instant when the message is printed out
      EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
      INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout,                                    ! Output of erhandler (ANSYS own subroutine for displaying program errors)                                 
     &  isubst,                                    ! Current substep number (sub-in)
     &  ieqitr,                                    ! Current equilibrium iteration number (in)
     &  elId                                       ! Element number (sub-in)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!     Calculate the crack growth rate based on the energy release rate through Paris' law
      ! Obtain the total energy release rate and mode mixity from the J-integral calculations
      G = GI + GII
      modemixity = GII/G
      ! In order to avoid errors due to a spurious very small negative GII, which would make the mode mixity equal to -0., let's make it zero
      ! This is done because later the mode mixity is a fractional power in the calculation of BKeta and a negative non-integer value would break it
      IF (modemixity < ZERO) THEN
            modemixity = ZERO
      END IF
      ! Calculate the parameters of the Paris law for the current mode mixity. Equation (44) from [2], removing the logs from teh calculation of C     
      cParis = ((cParisII/(cParisI*cParism))**((modemixity)**2))* 
     & (cParism**(modemixity))* cParisI
      mParis = (modemixity)**2*(mParisII-mParisI-mParism) 
     & + modemixity*mParism + mParisI
      ! Obtain the global critical energy release rate through the BK criterion. Equation (45) from [2]
      GcGlobal = GIc + (GIIc-GIc)*(modemixity)**BKeta
      ! Use the Paris' law to obtain the crack growth rate. Equation (43) from [2], but without using the threshold Gth and using Gmax instead of DeltaG
      dadN = cParis*(G/GcGlobal)**mParis     
!
!     Using the target crack growth per substep, calculate the cycles increment for the next crack propagation
      DeltaNProp = Deltaatarget/dadN 
!
! ! Debugging
!       errout = erhandler('USERELEM',1003,2,
!      & 'Debugging CyclePropagation.
!      & ERR: %F, %F, %F, %F.  
!      & Calculations: %F, %F, %F, %F, %F, %F, %F.  
!      & DeltaN: %F, %F.',
!      & (/GI, GII, G, modemixity,
!      & cParis, mParis, GcGlobal, dadN,
!      & Deltaatarget, GIc, GIIc, 
!      & DeltaNProp, DeltaNProp/1E6 /), ' ')
!
      RETURN
      END SUBROUTINE CyclesPropagation
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE AdjustCyclesPropagation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     Subroutine to adjust the number of cycles for propagation based on the number of iterations to aid convergence
      SUBROUTINE AdjustCyclesPropagation(DeltaNPropCrack, 
     & DeltaN, isubst, ieqitr, elId)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
        DOUBLE PRECISION
     &  DeltaNPropCrack,                           ! Increment of cycles based on the crack propagation calculated from the crack growth rate and the crack growth target per substep. Minimum of all cracks (sub-in)
     &  DeltaNtarget,                              ! Increment of cycles to obtain the target crack growth (inside) 
     &  DeltaN                                     ! Increment of cycles used in the calculations. Adjusted to aid convergence based on number of iterations (sub-out)
!    Next, the variables used in debugging are declared. Including some variables used to identify the instant when the message is printed out
      EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
      INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout,                                    ! Output of erhandler (ANSYS own subroutine for displaying program errors)                                 
     &  isubst,                                    ! Current substep number (sub-in)
     &  ieqitr,                                    ! Current equilibrium iteration number (sub-in)
     &  elId                                       ! Element number (sub-in)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!     Obtain the number of cycles for the target crack growth increment. Equation (28) from [2]
      DeltaNtarget = DeltaNPropCrack
!
!     Adjust the cycle increment to achieve the required convergence 
      ! Decrease slightly the number of cycles for iterations in between 30 and 39
      IF (INT(ieqitr/10) > 2 .AND.  
     & INT(ieqitr/10) < 4) THEN
            DeltaN = 0.99*DeltaNtarget
      ! Increase slightly the number of cycles for iterations in between 40 and 59
      ELSEIF (INT(ieqitr/10) > 3 .AND. 
     & INT(ieqitr/10) < 6) THEN
            DeltaN = 1.01*DeltaNtarget
      ! Decrease a bit more the number of cycles for iterations in between 60 and 79
      ELSEIF (INT(ieqitr/10) > 5 .AND. 
     & INT(ieqitr/10) < 8) THEN
            DeltaN = 0.98*DeltaNtarget
      ! Increase a bit more the number of cycles for iterations in between 80 and 99
      ELSEIF (INT(ieqitr/10) > 7 .AND. 
     & INT(ieqitr/10) < 10) THEN
            DeltaN = 1.02*DeltaNtarget
      ! Decrease a bit more the number of cycles for iterations in between 100 and 119
      ELSEIF (INT(ieqitr/10) > 9 .AND. 
     & INT(ieqitr/10) < 12) THEN
            DeltaN = 0.97*DeltaNtarget
      ! Increase a bit more the number of cycles for iterations in between 120 and 139
      ELSEIF (INT(ieqitr/10) > 11 .AND. 
     & INT(ieqitr/10) < 14) THEN
            DeltaN = 1.03*DeltaNtarget
      ! Decrease a bit more the number of cycles for iterations in between 140 and 159
      ELSEIF (INT(ieqitr/10) > 13 .AND. 
     & INT(ieqitr/10) < 16) THEN
            DeltaN = 0.95*DeltaNtarget
      ! Increase a bit more the number of cycles for iterations in between 160 and 179
      ELSEIF (INT(ieqitr/10) > 15 .AND. 
     & INT(ieqitr/10) < 18) THEN
            DeltaN = 1.05*DeltaNtarget
      ! Decrease a bit more the number of cycles for iterations in between 180 and 199
      ELSEIF (INT(ieqitr/10) > 17 .AND. 
     & INT(ieqitr/10) < 20) THEN
            DeltaN = 0.9*DeltaNtarget
      ! Increase a bit more the number of cycles for iterations in between 200 and 219
      ELSEIF (INT(ieqitr/10) > 19 .AND. 
     & INT(ieqitr/10) < 22) THEN
            DeltaN = 1.1*DeltaNtarget
      ! Decrease a bit more the number of cycles for iterations in between 220 and 259
      ELSEIF (INT(ieqitr/10) > 21 .AND. 
     & INT(ieqitr/10) < 26) THEN
            DeltaN = 0.8*DeltaNtarget
      ! Increase a bit more the number of cycles for iterations in between 260 and 309
      ELSEIF (INT(ieqitr/10) > 25 .AND. 
     & INT(ieqitr/10) < 31) THEN
            DeltaN = 1.2*DeltaNtarget
      ! Decrease a bit more the number of cycles for iterations in between 310 and 349
      ELSEIF (INT(ieqitr/10) > 30 .AND. 
     & INT(ieqitr/10) < 35) THEN
            DeltaN = 0.7*DeltaNtarget
      ! Increase a bit more the number of cycles for iterations in between 350 and 409
      ELSEIF (INT(ieqitr/10) > 34 .AND. 
     & INT(ieqitr/10) < 41) THEN
            DeltaN = 1.3*DeltaNtarget
      ! Decrease a bit more the number of cycles for iterations in between 410 and 449
      ELSEIF (INT(ieqitr/10) > 40 .AND. 
     & INT(ieqitr/10) < 45) THEN
            DeltaN = 0.5*DeltaNtarget
      ! Increase a bit more the number of cycles for iterations in between 439 and 509
      ELSEIF (INT(ieqitr/10) > 44 .AND. 
     & INT(ieqitr/10) < 51) THEN
            DeltaN = 1.5*DeltaNtarget
      ! Increase a bit more the number of cycles for iterations greater than 510
      ELSEIF (INT(ieqitr/10) > 50) THEN
            DeltaN = 0.01*DeltaNtarget
      ELSE 
            DeltaN = DeltaNtarget
      END IF 
!
! !     Debugging
!       IF (isubst == 10) THEN
!             errout = erhandler('USERELEM',1004,2,
!      & 'Debugging Modification DeltaN
!      & Substep %F, iteration %F and element %F
!      & Values: %F, %F, %F. ',
!      & (/isubst*ONE, ieqitr*ONE, elId*ONE,
!      & DeltaN, DeltaNtarget, ieqitr*ONE/), ' ')
!       END IF
!
!     End of the subroutine
      RETURN
      END SUBROUTINE AdjustCyclesPropagation
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE FatigueInitiation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     Subroutine to determine whether initiation takes place
      ! This subroutine first identifies the elements in five possible 'initiation cracks' based on the elements most likely to initiate 
      ! The minimum cycle increment for the five provitional 'initiation cracks' is compared to the one from propagation, making DeltaN equal to the lowest of the two
      ! If initiation takes place, the elements belonging to the 'initiation crack' are trhen flagged
      SUBROUTINE FatigueInitiation(InitiationLength, ElemLengths, 
     & Neighbours, InitiatedCrack, DeltaNInitElem, DeltaNInitCrack,
     & SeverityInitElem, SlopeInitElem, DInitElem, DeltaN,  
     & isubst, ieqitr, elId)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      DOUBLE PRECISION
     &  SumLengths,                                ! Sum of the lengths of the elements in the initiated crack (inside)
     &  InitiationLength,                          ! Length of the crack to be introduced when initiation is detected (sub-in)
     &  ElemLengths(1e5),                          ! Variable containing the element length for every element (sub-in)
     &  DeltaNInitElem(1e5),                       ! Increment of cycles based on the crack initiation for each elId (sub-in)
     &  CopyDeltaNInitElem(1e5),                   ! Copy of the previous variable, which will be modified at the end of each iteration to identify succesive minimum values (inside)
     &  SeverityInitElem(1e5),                     ! Variable containing the severity, ratio of current traction over static strength, of each element (sub-in)
     &  SlopeInitElem(1e5),                        ! Variable containing the initiation slope of each element based on their displacement-based mode mixity (sub-in)
     &  DInitElem(1e5),                            ! Variable containing the initiation damage variable of each element (sub-in)
     &  SumSeverities,                             ! Sum of the severities of the initiated elements for the calculation of DeltaNInitCrack (inside)
     &  SumSlopes,                                 ! Sum of the slopes of the initiated elements for the calculation of DeltaNInitCrack (inside)
     &  SumDis,                                    ! Sum of the initiation damage variables of the initiated elements for the calculation of DeltaNInitCrack (inside)
     &  ProvDeltaNInitCrack(5),                    ! Increment of cycles to initiation for the five provisional initiated cracks (inside)
     &  DeltaNInitCrack,                           ! Increment of cycles to initiation for the whole crack (sub-out)
     &  DeltaN                                     ! Increment of cycles used in the calculations (sub-inout)
      INTEGER
     &  InitiatedElem(1),                          ! Element with the lowest number of cycles to initiation. MINLOC always returns an array (inside)  
     &  ProvInitiatedCrack(5, 1e5),                ! Variable indicating whether this element belongs provisionally to an initiated crack (inside) 
     &  InitiatedCrack(1e5),                       ! Variable indicating whether this element belongs to an initiated crack (sub-out) 
     &  NeighbourElem,                             ! Neighbour element number. Used when identifying the elements belonging to each crack (inside)
     &  Neighbours(1e5, 4),                        ! Variable containing the nodes to left (i,1) and right (i,2) of element i, and the damage status in the last integration point (i,3) (sub-in)
     &  LastHigherNeighbour,                       ! Last statically damaged element at a higher coordinate than InitiatedElem (inside)
     &  LastLowerNeighbour,                        ! Last statically damaged element at a higher coordinate than InitiatedElem (inside)
     &  ToHigher,                                  ! Variable indicating if 1 that elements at a higher coordinate should be added because there are no more neighbours on the other side (inside)
     &  ToLower,                                   ! Variable indicating if 1 that elements at a lower coordinate should be added because there are no more neighbours on the other side (inside)
     &  NumInitCrack(1),                           ! Number of the provisional initiated crack with teh lowest cycle increment. Different from the other crack numbering (inside)
     &  i,                                         ! Used in loops (inside)
     &  j,                                         ! Used in loops (inside)
     &  k                                          ! Used in loops (inside)
!    Next, the variables used in debugging are declared
      EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
      INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout,                                    ! Output of erhandler (ANSYS own subroutine for displaying program errors)                                 
     &  isubst,                                    ! Current substep number (sub-in)
     &  ieqitr,                                    ! Current equilibrium iteration number (in)
     &  elId                                       ! Element number (sub-in)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!     Zero the flag indicating whether the elements belong provisionally to the initiated cracks
      ProvInitiatedCrack(:, :) = ZERO
      ! Create a copy of the vector with all the initiation cycle increments. 
      ! This variable will be modified to identify succesive minimum values by setting an artificially high value in the previous minimum
      CopyDeltaNInitElem = DeltaNInitElem
!
!     Loop to find the provisional initiated cracks from the five elements most likely to go through initiation next
      MainLoop: DO i = 1, 5
!           Find the element with the lowest number of cycles to initiation
            InitiatedElem = MINLOC(CopyDeltaNInitElem) 
            ! Exit the loop if the detected element has an artificially high DeltaNInit or it is not a cohesive element (its length has not been calculated)
            IF (MINVAL(CopyDeltaNInitElem) 
     & > 1.0D10 - 1.0D-10 .OR. 
     & ElemLengths(InitiatedElem(1)) < 1.0D-10) THEN
                  EXIT MainLoop
            END IF
            ! Add this element to the provitional initiation crack
            ProvInitiatedCrack(i, InitiatedElem(1)) = 1
            ! Add its length to the sum of the initiation element lengths
            SumLengths = ElemLengths(InitiatedElem(1))
            !
!           Identify the neighbouring elements which are already statically damaged
            ! Loop of elements in the direction towards a higher x coordinate
            NeighbourElem = InitiatedElem(1)
            OverallLoopHigh: DO k = 1, 1e5
                  LoopOneElemHigh: DO j = 1, 1e5
                        ! Identify the neighbouring element based on coinciding nodes
                        IF (Neighbours(j, 1) == 
     & Neighbours(NeighbourElem, 2)) THEN
                              ! Only save this element and continue, if the element is already damaged (its stiffness damage variable Dk is different to 0 in one side or the other)
                              IF (Neighbours(j, 3) > 1.0D-10 
     & .OR. Neighbours(j, 4) > 1.0D-10) THEN
                                    ! Save this new element to start the next lap in the loop with it
                                    NeighbourElem = j
                                    ! Add this element to the provitional initiation crack, as well as its length to the sum of all lengths
                                    ProvInitiatedCrack(i,j) = 1
                                    SumLengths = 
     & SumLengths + ElemLengths(j)
                                    ! Exit and restart the loop from element j
                                    EXIT LoopOneElemHigh
                              ! Exit the overall loop if the next neighbouring element is not already damaged
                              ELSE
                                    EXIT OverallLoopHigh
                              END IF
                        END IF
                  END DO LoopOneElemHigh
                  ! Exit the loop if there are no more neighbouring elements
                  IF (j == 1e5+1) THEN
                        EXIT OverallLoopHigh
                  END IF
            END DO OverallLoopHigh
            ! Save the last element to have been identified as part of the 'initiation crack'
            LastHigherNeighbour = NeighbourElem
            ! Loop of elements in the direction towards a lower x coordinate
            NeighbourElem = InitiatedElem(1)
            OverallLoopLow: DO k = 1, 1e5
                  LoopOneElemLow: DO j = 1, 1e5
                        ! Identify the neighbouring element based on coinciding nodes
                        IF (Neighbours(j, 2) == 
     & Neighbours(NeighbourElem, 1)) THEN
                              ! Only save this element and continue, if the element is already damaged (its stiffness damage variable Dk is different to 0 in one side or the other)
                              IF (Neighbours(j, 3) > 1.0D-10 
     & .OR. Neighbours(j, 4) > 1.0D-10) THEN
                                    ! Save this new element to start the next lap in the loop with it
                                    NeighbourElem = j
                                    ! Add this element to the provitional initiation crack, as well as its length to the sum of all lengths
                                    ProvInitiatedCrack(i,j) = 1
                                    SumLengths = 
     & SumLengths + ElemLengths(j)
                                    EXIT LoopOneElemLow
                              ! Exit the overall loop if the next neighbouring element is not already damaged
                              ELSE
                                    EXIT OverallLoopLow
                              END IF
                        END IF
                  END DO LoopOneElemLow
                  ! Exit the loop if there are no more neighbouring elements
                  IF (j == 1e5+1) THEN
                        EXIT OverallLoopLow
                  END IF
            END DO OverallLoopLow
            ! Save the last element to have been identified as part of the 'initiation crack'
            LastLowerNeighbour = NeighbourElem
            !
            ! Next, add not already damaged elements to the crack, as long as the sum of the element lengths in the initiated crack is below the value in InitiationLength
            ! These elements are added based on where the number of cycles to initiation is lower 
            ! It also accounts for the possibility that there might be no more neighbouring elements to one side of the crack
            ! First, zero the variables, which if made 1 indicate that no more neighbouring elements can be found to one side of the crack
            ToHigher = 0
            ToLower = 0
            OverallLoop2: DO WHILE (SumLengths < InitiationLength)
                  ! Loop through the neighbouring elements to find those at higher x coordinates
                  ! If the element with the lowest DeltaNInit is to a higher coordinate and there are still neighbours on that side
                  ! Or if there are no more neighbouring elements at a lower coordinate (ToHigher = 1)
                  ! In the case that there are no statically damaged elements, LastLowerNeighbour and LastHigherNeighbour are InitiatedElem, so .LE. is used to compare the DeltaNInitElem
                  IF ((DeltaNInitElem(LastHigherNeighbour) 
     & .LE. DeltaNInitElem(LastLowerNeighbour)
     & .AND. ToLower == 0) .OR. ToHigher == 1) THEN
                        NeighbourElem = LastHigherNeighbour
                        LoopOneElemHigh2: DO j = 1, 1e5
                              ! Identify the neighbouring element based on coinciding nodes
                              IF (Neighbours(j, 1) == 
     & Neighbours(NeighbourElem, 2)) THEN
                                    ! Save this new element to start the next lap in the loop with it
                                    LastHigherNeighbour = j
                                    ! Add this element to the provitional initiation crack, as well as its length to the sum of all lengths
                                    ProvInitiatedCrack(i,j) = 1
                                    SumLengths = 
     & SumLengths + ElemLengths(j)
                                    EXIT LoopOneElemHigh2
                              END IF 
                        END DO LoopOneElemHigh2
                        ! If there are no more neighbouring elements, make ToLower equal to 1 so that the algorithm checks on the other side of the crack (lower coordinate)
                        IF (j == 1e5+1) THEN
                              ToLower = 1
                        END IF     
                  ! Loop through the neighbouring elements to find those at lower x coordinates
                  ! If the element with the lowest DeltaNInit is to a lower coordinate and there are still neighbours on that side
                  ! Or if there are no more neighbouring elements at a higher coordinate (ToLower = 1)
                  ! In the case that there are no statically damaged elements, LastLowerNeighbour and LastHigherNeighbour are InitiatedElem, so .LE. is used to compare the DeltaNInitElem
                  ELSE IF ((DeltaNInitElem(LastHigherNeighbour) 
     & > DeltaNInitElem(LastLowerNeighbour)
     & .AND. ToHigher == 0) .OR. ToLower == 1) THEN
                        NeighbourElem = LastLowerNeighbour
                        LoopOneElemLow2: DO j = 1, 1e5
                              ! Identify the neighbouring element based on coinciding nodes
                              IF (Neighbours(j, 2) == 
     & Neighbours(NeighbourElem, 1)) THEN
                                    ! Save this new element to start the next lap in the loop with it
                                    LastLowerNeighbour = j
                                    ! Add this element to the provitional initiation crack, as well as its length to the sum of all lengths
                                    ProvInitiatedCrack(i,j) = 1
                                    SumLengths = 
     & SumLengths + ElemLengths(j)
                                    EXIT LoopOneElemLow2
                              END IF 
                        END DO LoopOneElemLow2
                        ! If there are no more neighbouring elements, make ToHigher equal to 1 so that the algorithm checks on the other side of the crack (higher coordinate)
                        IF (j == 1e5+1) THEN
                              ToHigher = 1
                        END IF
                  END IF
            END DO OverallLoop2
            !
            ! Calculate for each provitional 'initiation crack' the increment cycle by averaging the elements' initiation slopes, severities and initiation damage variables
            ! First, zero the severities, slopes and initiation damage variables
            SumSeverities = ZERO
            SumSlopes = ZERO
            SumDis = ZERO
            ! Zero the variable k, which is used to sum the number of elements in each provitional 'initiation crack'
            k = 0
            ! Loop through all elements
            DeltaNLoop: DO j = 1, 1e5
                  ! Identify the elements provisionally in this 'initiation crack'
                  IF (ProvInitiatedCrack(i,j) == 1) THEN
                        ! Sum up the severities and initiation slopes
                        SumSeverities = SumSeverities +
     & SeverityInitElem(j) 
                        SumSlopes = SumSlopes +
     & SlopeInitElem(j) 
                        ! When summing the damage initiation variables, limit their value to 1 if they are greater
                        IF (DInitElem(j) > ONE) THEN
                              SumDis = SumDis + ONE
                        ELSE
                        SumDis = SumDis +
     & DInitElem(j) 
                        END IF
                        ! Use k to sum the amount of elements
                        k = k + 1
                  END IF 
            END DO DeltaNLoop
            ! Calculate the increment cycle for this 'initiation crack' by averaging the elements' initiation slopes, severities and initiation damage variables
            ProvDeltaNInitCrack(i) = (ONE-SumDis/k)* 
     & (10**((ONE-SumSeverities/k)/(SumSlopes/k)))
            !
            ! Give an extremely high value to the DeltaNInitElem for the current initiated element, so that in the next lap a new minimum will be found
            CopyDeltaNInitElem(InitiatedElem(1)) = 1.0D10
      END DO MainLoop
      !
      ! The final 'initiation crack' will be the one from all of the provisional ones with the lowest cycle increment
      DeltaNInitCrack = MINVAL(ProvDeltaNInitCrack)
      ! Zero the cycle increment in case it is negative. This can occur if this crack initiates after another and the load redistribution has increased greatly the severity for the second crack
      ! It is thus applied in this case that both cracks would originate at the same time, even if in different substeps in the analysis
      IF (DeltaNInitCrack < ZERO) THEN
            DeltaNInitCrack = ZERO
      END IF
      ! Save the number of the chosen provitional 'initiation crack'
      NumInitCrack = MINLOC(ProvDeltaNInitCrack)
      !
      ! If the calculated cycle increment according to initiation is lower than the one due to propagation
      IF (DeltaNInitCrack < DeltaN) THEN
            ! Initiation occurs and DeltaN is set to the initiation cycle increment
            DeltaN = DeltaNInitCrack
            DO j = 1, 1e5
                  ! Identify the elements provisionally in the 'initiation crack' and flag them with InitiatedCrack so that they will become fully damaged
                  IF (ProvInitiatedCrack(NumInitCrack(1),j) == 1) THEN
                        InitiatedCrack(j) = 1
                  END IF
            END DO
      END IF
      ! Debugging
      errout = erhandler('USERELEM',1005,2,
     & 'Debugging FatigueInitiation
     & Substep %F, iteration %F and element %F
     & Values: %F, %F.',
     & (/isubst*ONE, ieqitr*ONE, elId*ONE,
     & DeltaNInitCrack/1e6, NumInitCrack(1)*ONE
     & /), ' ')
!
!     End of the subroutine
      RETURN
      END SUBROUTINE FatigueInitiation
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE Slopes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     Subroutine to calculate the derivatives of the opening and shearing displacements wrt the crack growth (ddelta3da & ddelta1da)
      SUBROUTINE Slopes(xRef, nDim, nNodes, nlgeom, xCur, 
     & TotValDofs, nUsrDof, PropNode, ddelta3da, ddelta1da, 
     & isubst, ieqitr, elId)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      INTEGER ::
     &  nDim,                                      ! Number of dimensions of the problem (sub-in)
     &  nNodes,                                    ! Number of nodes of the element (sub-in)
     &  IntP,                                      ! Integration point number (inside)
     &  iInt,                                      ! Integration counter along natural coordinate eta (inside)
     &  jInt,                                      ! Integration counter along natural coordinate xi (inside)
     &  nlgeom,                                    ! Large deformation key [from nlgeom command] (sub-in)
     &  nUsrDof,                                   ! Number of DOFs of this element (matrix and load vector size) (sub-in)
     &  PropNode                                   ! Propagation node numbering: 0 indicates node 1 (lower  coordinate) and 1 the node 2 (greater x coordinate) (sub-in)
      DOUBLE PRECISION ::
     &  elLength,                                  ! Approximation of the element length (inside) 
     &  xCur(nDim,nNodes),                         ! Nodal coordinates in current configuration (sub-in)
     &  IntPCoords2(2**2),                         ! Natural coordinates for the integration points using Newton-Cotes and 2 integration points (inside)
     &  IntPWeights2(2**2),                        ! Weights for the integration points using Newton-Cotes and 2 integration points (inside)     
     &  ShapeFactors(nNodes),                      ! Shape function factors in local coordinates (inside)
     &  DerivShapeFactors(nDim-1,nNodes),          ! Factors for the derivatives of the shape functions (inside)
     &  xRef(nDim,nNodes),                         ! Nodal coordinates in initial configuration (sub-in)     
     &  Theta(nDim, nDim),                         ! Matrix of directional vectors. Orthogonal rotation tensor relating the local and fixed coordinate systems (inside)
     &  dA,                                        ! Area scale factor for current integration point value (inside)
     &  Bglobal(nDim, nUsrDof),                    ! B-matrix in global coordinates (inside)
     &  TotValDofs(nUsrDof),                       ! Total values of DOFs (displacements) from time 0 (sub-in)
     &  DisplJump2(2**2, nDim),                    ! Displacement jump at integration point using Newton-Cotes and 2 integration points (inside)
     &  PropDir,                                   ! Propagation direction: -1 for along +eta and 1 for along -eta (inside)
     &  ddelta3da,                                 ! Derivative of the opening displacement with respect to the crack growth a (sub-out)
     &  ddelta1da                                  ! Derivative of the mode II displacement with respect to the crack growth a (sub-out) 
!    Next, the variables used in debugging are declared. Including some variables used to identify the instant when the message is printed out
      EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
      INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout,                                    ! Output of erhandler (ANSYS own subroutine for displaying program errors)                                 
     &  isubst,                                    ! Current substep number (sub-in)
     &  ieqitr,                                    ! Current equilibrium iteration number (in)
     &  elId                                       ! Element number (sub-in)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
! First, approximate the element length as the distance at the mid-surface between nodes 1 and 2
      elLength = SQRT( ((xCur(1, 2)+xCur(1, 6))/2 
     & - (xCur(1, 1)+xCur(1, 5))/2)**2 
     & +  ((xCur(2, 2)+xCur(2, 6))/2 
     & - (xCur(2, 1)+xCur(2, 5))/2)**2)
!
! Next, the displacement jump must be recalculated at both extremes of the element's length, using thus 2 integration points and Newton-Cotes so that the integration points coincide with the nodes
      CALL IntegrationPoints(2, IntPCoords2, IntPWeights2, 1)
      IntP = 0 
      DO iInt = 1, 2
            DO jInt = 1, 2
                  IntP = IntP + 1
                  CALL ShapeFunc(IntPCoords2(iInt), IntPCoords2(jInt), 
     &  ShapeFactors, DerivShapeFactors, nNodes, nDim)
                  CALL Jacobian(nDim, nNodes, DerivShapeFactors, xRef, 
     &  nlgeom, xCur, Theta, dA, 2, IntPWeights2, iInt, jInt)
                  CALL Displacements(nDim, nUsrDof, nNodes,   
     &  ShapeFactors, Bglobal, Theta, DisplJump2, 2, IntP, TotValDofs)
            END DO
      END DO
!
! Once the displacement jumps have been recalculated, the derivatives can be obtained
      ! The direction of propagation is imposed
      IF (PropNode == 0) THEN
            PropDir = ONE
      ELSE
            PropDir = -ONE 
      END IF
      ! The derivatives of the displacement jumps with respect to the crack growth a are approximated under the assumption of self-similar growth to the derivatives of the displacement jumps with respect to
      ! the crack growth coordinate xa, but with opposite signs. This is a correction of equations (25) and (26) from [2].
      ! The derivatives of the displacement jumps with respect to the crack growth coordinate xa can then be calculated as the difference of the displacement jumps at both extremes of the element's length 
      ! over such length
      ! Equation (38) in [2], which is incorrect if used to obtain the derivative with respect to xa as indicated in the paper, but it can be used to obtain the deirvative wrt a
      ddelta3da = PropDir*(DisplJump2(3,3) - 
     & DisplJump2(1,3)) / elLength
      ddelta1da = PropDir*(DisplJump2(3,1) - DisplJump2(1,1))
     & / elLength
!
! Debugging
!       errout = erhandler('USERELEM',1006,2,
!      & 'ERROR in ddelta3da %F
!      & Substep %F, iteration %F and element %F
!      & PropDir: %F.
!      & In detail: %F, %F.
!      & Others: %F, %F, %F, %F. ',
!      & (/ddelta3da,
!      & isubst*ONE, ieqitr*ONE, elId*ONE,
!      & PropDir, 
!      & DisplJump2(3,3), DisplJump2(1,3),
!      & elLength/), ' ')
!
      RETURN
      END SUBROUTINE Slopes
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE IntegrationPoints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     Subroutine to obtain the coordinates and weights of the integration points
      SUBROUTINE IntegrationPoints(NIntP, IntPCoords, IntPWeights, 
     &  IntMethodFlag)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      INTEGER 
     &  NIntP,                                     ! Number of integration points (sub-in)
     &  IntMethodFlag,                             ! Flag to select the integration method: Newton-Cotes (1) or Gauss (2) (sub-in)
     &  i                                          ! Used in loops (inside)
      DOUBLE PRECISION
     &  IntPCoords(NIntP),                         ! Natural coordinates for the integration points (sub-out)
     &  IntPWeights(NIntP),                        ! Weights for the integration points (sub-out)
     &  hNewton                                    ! Distance between integration points in Newton-Cotes integration (inside)
!    Next, the variables used in debugging are declared
        EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
        INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout                                     ! Output of erhandler (ANSYS own subroutine for displaying program errors) (inside)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!     Select the integration method
      !
      ! Newton-Cotes integration
      IF (IntMethodFlag == 1) THEN
         ! h is calculated as (b-a)/NIntP where b=1 and a=-1 in this case
         hNewton = TWO/((NIntP-1)*ONE)
         ! The weights and coordinates are then defined
         DO i = 1, NIntP
            ! For the extremes, the weights differ
            IF (i == 1 .OR. i == NIntP) THEN
               IntPWeights(i) = HALF*hNewton
            ELSE 
               IntPWeights(i) = hNewton
            END IF
            ! The coordinates vary from -1 to 1, spaced by h
            IntPCoords(i) = -ONE + hNewton*(i-1)
         END DO
      !
      ! Gauss integration
      ELSEIF (IntMethodFlag == 2) THEN
          ! Each case represents a certain number of Gauss integration points
          SELECT CASE (NIntP)
            CASE(2)
                  IntPCoords(1:2) = (/-.5773502693, .5773502693/)
                  IntPWeights(1:2) = (/1.00000000, 1.00000000/)
            CASE(3)
                  IntPCoords(1:3) = (/-0.7745966692,0.,0.7745966692/)
                  IntPWeights(1:3) = (/.5555555554, .8888888888, 
     & 0.5555555554/)
            CASE(4)
                  IntPCoords(1:4) = (/-0.8611363114,-0.3399810437,
     & 0.3399810437,0.8611363114/)
                  IntPWeights(1:4) = (/0.3478548434,0.6521451552,
     & 0.6521451552,0.3478548434/)
            CASE(5)
                  IntPCoords(1:5) = (/-0.9061798457,-0.5384693100,0.,
     & 0.5384693100,0.9061798457/)
                  IntPWeights(1:5) = (/0.2369268822,0.4786286702,
     & 0.5688888890,0.4786286702,0.2369268822/)
            CASE(6)
                  IntPCoords(1:6) = (/-0.9324695139,-0.6612093866,
     & -0.2386191855,0.2386191855,0.6612093866,0.9324695139/)
                  IntPWeights(1:6) = (/0.1713244879,0.3607615738,
     & 0.4679139338,0.4679139338,0.3607615738,0.1713244879/)
            CASE(7)
                  IntPCoords(1:7) = (/-0.9491079124,-0.7415311855,
     & -0.4058451510,0.,0.4058451510,0.7415311855,0.9491079124/)
                  IntPWeights(1:7) = (/0.1294849672,0.2797053910,
     & 0.3818300496,0.4179591836,0.3818300496,0.2797053910,
     & 0.1294849672/)
            CASE(8)
                  IntPCoords(1:8) = (/-0.9602898564,-0.7966664774,
     & -0.5255324098,-0.1834346428,0.1834346428,0.5255324098,
     & 0.7966664774,0.9602898564/)
                  IntPWeights(1:8) = (/0.1012285344,0.2223810342,
     & 0.3137066454,0.3626837836,0.3626837836,0.3137066454,
     & 0.2223810342,0.1012285344/)
            CASE(9)
                  IntPCoords(1:9) = (/-0.9681602395,-0.8360311074,
     & -0.6133714327,-0.3242534235,0.,0.3242534235,0.6133714327,
     & 0.8360311074,0.9681602395/)
            IntPWeights(1:9) = (/0.08127438834,0.1806481615,
     & 0.2606106964,0.3123470772,0.3302393550,0.3123470772,
     & 0.2606106964,0.1806481615,0.08127438834/)
            CASE(10)
                  IntPCoords(1:10) = (/-0.9739065285,-0.8650633667,
     & -0.6794095683,-0.4333953942,-0.1488743390,0.1488743390,
     & 0.4333953942,0.6794095683,0.8650633667,0.9739065285/)
            IntPWeights(1:10) = (/0.06667134380,0.1494513492,
     & 0.2190863626,0.2692667196,0.2955242248,0.2955242248,
     & 0.2692667196,0.2190863626,0.1494513492,0.06667134380/)
            CASE(15)
                  IntPCoords(1:15) = (/-0.9879925180,-0.9372733924,
     & -0.8482065834,-0.7244177314,-0.5709721726,-0.3941513471,
     & -0.2011940940,0.,0.2011940940,0.3941513471,0.5709721726,
     & 0.7244177314,0.8482065834,0.9372733924,0.9879925180/)
                  IntPWeights(1:15) = (/0.03075324124,0.07036604748,
     & 0.1071592204,0.1395706781,0.1662692057,0.1861610000,
     & 0.1984314853,0.2025782420,0.1984314853,0.1861610000,
     & 0.1662692057,0.1395706781,0.1071592204,0.07036604748,
     & 0.03075324124/)
            CASE(20)
                  IntPCoords(1:20) = (/-0.9931285992,-0.9639719273,
     & -0.9122344283,-0.8391169718,-0.7463319065,-0.6360536808,
     & -0.5108670020,-0.3737060887,-0.2277858511,-0.07652652113,
     & 0.07652652113,0.2277858511,0.3737060887,0.5108670020,
     & 0.6360536808,0.7463319065,0.8391169718,0.9122344283,
     & 0.9639719273,0.9931285992/)
                  IntPWeights(1:20) = (/0.01761400797,0.04060143030,
     & 0.06267204900,0.08327674140,0.1019301200,0.1181945322,
     & 0.1316886387,0.1420961093,0.1491729865,0.1527533872,
     & 0.1527533872,0.1491729865,0.1420961093,0.1316886387,
     & 0.1181945322,0.1019301200,0.08327674140,0.06267204900,
     & 0.04060143030,0.01761400797/)
            CASE(30)
                  IntPCoords(1:30) = (/-0.9968934841,-0.9836681233,
     & -0.9600218650,-0.9262000475,-0.8825605358,-0.8295657624,
     & -0.7677774321,-0.6978504948,-0.6205261830,-0.5366241482,
     & -0.4470337695,-0.3527047255,-0.2546369262,-0.1538699136,
     & -0.05147184255,0.05147184255,0.1538699136,0.2546369262,
     & 0.3527047255,0.4470337695,0.5366241482,0.6205261830,
     & 0.6978504948,0.7677774321,0.8295657624,0.8825605358,
     & 0.9262000475,0.9600218650,0.9836681233,0.9968934841/)
                  IntPWeights(1:30) = (/0.007968194496,0.01846646900,
     & 0.02878470854,0.03879919362,0.04840267292,0.05749315636,
     & 0.06597422980,0.07375597474,0.08075589534,0.08689978736,
     & 0.09212252218,0.09636873720,0.09959342068,0.1017623896,
     & 0.1028526529,0.1028526529,0.1017623896,0.09959342068,
     & 0.09636873720,0.09212252218,0.08689978736,0.08075589534,
     & 0.07375597474,0.06597422980,0.05749315636,0.04840267292,
     & 0.03879919362,0.02878470854,0.01846646900,0.007968194496/)
            CASE(50)
                  IntPCoords(1:50) = (/-0.9988664044,-0.9940319695,
     & -0.9853540841,-0.9728643851,-0.9566109552,-0.9366566189,
     & -0.9130785566,-0.8859679795,-0.8554297694,-0.8215820709,
     & -0.7845558329,-0.7444943022,-0.701552468706823,-0.6558964657,
     & -0.607702927184950,-0.557158304514650,-0.504458144907465,
     & -0.449806334974039,-0.393414311897565,-0.335500245419437,
     & -0.276288193779532,-0.216007236876042,-0.154890589998146,
     & -0.0931747015600861,-0.0310983383271889,0.0310983383271889,
     & 0.0931747015600861,0.154890589998146,0.216007236876042,
     & 0.276288193779532,0.335500245419437,0.393414311897565,
     & 0.449806334974039,0.504458144907465,0.557158304514650,
     & 0.607702927184950,0.6558964657,0.701552468706823,
     & 0.7444943022,0.7845558329,0.8215820709,0.8554297694,
     & 0.8859679795,0.9130785566,0.9366566189,0.9566109552,
     & 0.9728643851,0.9853540841,0.9940319695,0.9988664044/)
                  IntPWeights(1:50) = (/0.002908620030,0.006759803022,
     & 0.01059055023,0.01438082259,0.01811555984,0.02178024244,
     & 0.02536067280,0.02884299330,0.03221372792,0.03545983598,
     & 0.03856875662,0.04152846292,0.0443275043407263,0.04695505138,
     & 0.0494009384516935,0.0516557030373905,0.0537106218784147,
     & 0.0555577448137462,0.0571899256357255,0.0586008498273267,
     & 0.0597850587471728,0.0607379708565016,0.0614558995460401,
     & 0.0619360674583825,0.0621766166808138,0.0621766166808138,
     & 0.0619360674583825,0.0614558995460401,0.0607379708565016,
     & 0.0597850587471728,0.0586008498273267,0.0571899256357255,
     & 0.0555577448137462,0.0537106218784147,0.0516557030373905,
     & 0.0494009384516935,0.04695505138,0.0443275043407263,
     & 0.04152846292,0.03856875662,0.03545983598,0.03221372792,
     & 0.02884299330,0.02536067280,0.02178024244,0.01811555984,
     & 0.01438082259,0.01059055023,0.006759803022,0.002908620030/)
	      CASE(70)					
                  IntPCoords(1:70) = (/-0.9994182860,-0.9969362520,
     & -0.9924760552,-0.9860455581,-0.9776574059,-0.9673282237,
     & -0.9550785091,-0.9409325790,-0.9249185169,-0.9070681162,
     & -0.8874168169,-0.8660036342,-0.8428710820,-0.8180650876,
     & -0.7916349010,-0.7636329968,-0.7341149701,-0.7031394261,
     & -0.6707678641,-0.6370645546,-0.6020964125,-0.5659328637,
     & -0.5286457076,-0.4903089745,-0.4509987784,-0.4107931660,
     & -0.3697719616,-0.3280166093,-0.2856100105,-0.2426363595,
     & -0.1991809764,-0.1553301379,-0.1111709058,-0.06679095417,
     & -0.02227839529,0.02227839529,0.06679095417,0.1111709058,
     & 0.1553301379,0.1991809764,0.2426363595,0.2856100105,
     & 0.3280166093,0.3697719616,0.4107931660,0.4509987784,
     & 0.4903089745,0.5286457076,0.5659328637,0.6020964125,
     & 0.6370645546,0.6707678641,0.7031394261,0.7341149701,
     & 0.7636329968,0.7916349010,0.8180650876,0.8428710820,
     & 0.8660036342,0.8874168169,0.9070681162,0.9249185169,
     & 0.9409325790,0.9550785091,0.9673282237,0.9776574059,
     & 0.9860455581,0.9924760552,0.9969362520,0.9994182860/)
                  IntPWeights(1:70) = (/0.001492726023,0.003471897902,
     & 0.005447111300,0.007411770470,0.009361761102,0.01129318548,
     & 0.01320219054,0.01508498780,0.01693783641,0.01875705628,
     & 0.02053903828,0.02228024030,0.02397720790,0.02562657068,
     & 0.02722505476,0.02876948616,0.03025679824,0.03168403764,
     & 0.03304837226,0.03434709198,0.03557761898,0.03673750960,
     & 0.03782446126,0.03883631628,0.03977106556,0.04062685302,
     & 0.04140197902,0.04209490554,0.04270425672,0.04322882256,
     & 0.04366756144,0.04401960242,0.04428424654,0.04446096842,
     & 0.04454941716,0.04454941716,0.04446096842,0.04428424654,
     & 0.04401960242,0.04366756144,0.04322882256,0.04270425672,
     & 0.04209490554,0.04140197902,0.04062685302,0.03977106556,
     & 0.03883631628,0.03782446126,0.03673750960,0.03557761898,
     & 0.03434709198,0.03304837226,0.03168403764,0.03025679824,
     & 0.02876948616,0.02722505476,0.02562657068,0.02397720790,
     & 0.02228024030,0.02053903828,0.01875705628,0.01693783641,
     & 0.01508498780,0.01320219054,0.01129318548,0.009361761102,
     & 0.007411770470,0.005447111300,0.003471897902,0.001492726023/)
	      CASE(90)					
                  IntPCoords(1:90) = (/-0.9996469713,-0.9981403800,
     & -0.9954318121,-0.9915239288,-0.9864213650,-0.9801302513,
     & -0.9726581621,-0.9640140982,-0.9542084739,-0.9432531036,
     & -0.9311611875,-0.9179472951,-0.9036273479,-0.8882186005,
     & -0.8717396189,-0.8542102591,-0.8356516426,-0.8160861309,
     & -0.7955372992,-0.7740299070,-0.7515898690,-0.7282442239,
     & -0.7040211012,-0.6789496880,-0.6530601932,-0.6263838118,
     & -0.5989526868,-0.5707998704,-0.5419592845,-0.5124656801,
     & -0.4823545944,-0.4516623089,-0.4204258056,-0.3886827219,
     & -0.3564713059,-0.3238303696,-0.2907992431,-0.2574177260,
     & -0.2237260407,-0.1897647829,-0.1555748733,-0.1211975081,
     & -0.08667410942,-0.05204627514,-0.01735572915,0.01735572915,
     & 0.05204627514,0.08667410942,0.1211975081,0.1555748733,
     & 0.1897647829,0.2237260407,0.2574177260,0.2907992431,
     & 0.3238303696,0.3564713059,0.3886827219,0.4204258056,
     & 0.4516623089,0.4823545944,0.5124656801,0.5419592845,
     & 0.5707998704,0.5989526868,0.6263838118,0.6530601932,
     & 0.6789496880,0.7040211012,0.7282442239,0.7515898690,
     & 0.7740299070,0.7955372992,0.8160861309,0.8356516426,
     & 0.8542102591,0.8717396189,0.8882186005,0.9036273479,
     & 0.9179472951,0.9311611875,0.9432531036,0.9542084739,
     & 0.9640140982,0.9726581621,0.9801302513,0.9864213650,
     & 0.9915239288,0.9954318121,0.9981403800,0.9996469713/)
                  IntPWeights(1:90) = (/0.0009059354934,0.002107785028,
     & 0.003308869966,0.004506123082,0.005697979396,0.006882981816,
     & 0.008059695198,0.009226697606,0.01038258268,0.01152595708,
     & 0.01265544583,0.01376968559,0.01486733268,0.01594706794,
     & 0.01700758669,0.01804761297,0.01906589367,0.02006120030,
     & 0.02103233622,0.02197812928,0.02289743976,0.02378916154,
     & 0.02465221882,0.02548557252,0.02628821738,0.02705918730,
     & 0.02779755346,0.02850242536,0.02917295354,0.02980833148,
     & 0.03040779240,0.03097061398,0.03149611872,0.03198367296,
     & 0.03243268958,0.03284262702,0.03321299198,0.03354333758,
     & 0.03383326626,0.03408242838,0.03429052386,0.03445730190,
     & 0.03458256166,0.03466615210,0.03470797248,0.03470797248,
     & 0.03466615210,0.03458256166,0.03445730190,0.03429052386,
     & 0.03408242838,0.03383326626,0.03354333758,0.03321299198,
     & 0.03284262702,0.03243268958,0.03198367296,0.03149611872,
     & 0.03097061398,0.03040779240,0.02980833148,0.02917295354,
     & 0.02850242536,0.02779755346,0.02705918730,0.02628821738,
     & 0.02548557252,0.02465221882,0.02378916154,0.02289743976,
     & 0.02197812928,0.02103233622,0.02006120030,0.01906589367,
     & 0.01804761297,0.01700758669,0.01594706794,0.01486733268,
     & 0.01376968559,0.01265544583,0.01152595708,0.01038258268,
     & 0.009226697606,0.008059695198,0.006882981816,0.005697979396,
     & 0.004506123082,0.003308869966,0.002107785028,0.0009059354934/) 
            CASE DEFAULT  
                  ! Generate an error message if the number of integration points is not in one of the previous cases
                  errout = erhandler('USERELEM',5001,3,
     & 'ERROR. When using Gauss integration, the number of
     & integration points must be 2, 3, 4, 5, 6, 7, 8, 9, 10, 15,
     & 20, 30, 50, 70 or 90.', , ' ')   
          END SELECT
      ! Error if the flag has any other value
      ! ELSE
      END IF
!
! !     Debugging
!       errout = erhandler('USERELEM',1007,2,
!      & 'IntegrationPoints. Printing some values   
!      & %F, %F, %F, %F, %F.',
!      & (/IntMethodFlag*ONE, IntPWeights(1), IntPCoords(1),  
!      & IntPWeights(2), IntPCoords(2)/), ' ')
!
!     End of the subroutine
      RETURN
      END SUBROUTINE IntegrationPoints
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE ShapeFunc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     Subroutine to evaluate the shape functions and their derivatives for one integration point
      SUBROUTINE ShapeFunc(eta, xi, ShapeFactors, DerivShapeFactors,
     &  nNodes, nDim)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      INTEGER 
     &  nDim,                                      ! Number of dimensions of the problem (sub-in)
     &  nNodes,                                    ! Number of nodes of the element (sub-in)     
     &  i,                                         ! Used in loops (inside)
     &  j                                          ! Used in loops (inside)
        DOUBLE PRECISION
     &  eta,                                       ! Natural coordinate eta (x if interface along the x coordinate) (sub-in)
     &  xi,                                        ! Natural coordinate xi (y if the y axis is normal to the interface) (sub-in)
     &  ShapeFactors(nNodes),                      ! Shape function factors in local coordinates (sub-out) 
     &  DerivShapeFactors(nDim-1,nNodes)           ! Factors for the derivatives of the shape functions (sub-out)
!    Next, the variables used in debugging are declared
        EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
        INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout                                     ! Output of erhandler (ANSYS own subroutine for displaying program errors) (inside)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!     Initialise the variables
      DO i = 1, nNodes 
         ShapeFactors(i) = ZERO 
         DO j = 1, nDim-1
            DerivShapeFactors(j,i) = ZERO
         END DO
      END DO
!
!     Evaluate the shape function and their derivatives factors for an 8-node element in this integration point
      IF (nNodes .EQ. 8) THEN
         ! First four shape function factors. Using standard Lagrangian functions
         ShapeFactors(1) = (ONE-eta)*(ONE-xi)/FOUR     ! N1
         ShapeFactors(2) = (ONE+eta)*(ONE-xi)/FOUR     ! N2
         ShapeFactors(3) = (ONE+eta)*(ONE+xi)/FOUR     ! N3
         ShapeFactors(4) = (ONE-eta)*(ONE+xi)/FOUR     ! N4
         ! First eight shape function derivative factors
         ! DerivShapeFactors(1,i) contains the derivatives wrt eta and (2,i) wrt xi
         DerivShapeFactors(1,1) = -(ONE-xi)/FOUR     ! dN1/deta
         DerivShapeFactors(2,1) = -(ONE-eta)/FOUR    ! dN1/dxi
         DerivShapeFactors(1,2) = (ONE-xi)/FOUR      ! dN2/deta
         DerivShapeFactors(2,2) = -(ONE+eta)/FOUR    ! dN2/dxi
         DerivShapeFactors(1,3) = (ONE+xi)/FOUR      ! dN3/deta 
         DerivShapeFactors(2,3) = (ONE+eta)/FOUR     ! dN3/dxi
         DerivShapeFactors(1,4) = -(ONE+xi)/FOUR     ! dN4/deta
         DerivShapeFactors(2,4) = (ONE-eta)/FOUR     ! dN4/dxi
      ! Generate an error message if the number of elements was not set to 8
      ELSE
         errout = erhandler('USERELEM',5002,3,
     & 'ERROR. The number of nodes must be set to 8.', , ' ')
      END IF
      ! Modify the shape function factors and calculate the last eight shape function derivative factors
      DO i = 1, nNodes/2
         ! Last four shape function factors (top nodes) are the same as the ones previously defined
         ShapeFactors(i+nNodes/2) = ShapeFactors(i)
         ! First four shape functions (bottom nodes) are defined then as the negative versions 
         ! Expression (60) in [1] or original expression (4) in [4]
         ShapeFactors(i) = -ShapeFactors(i+nNodes/2)
         DO j = 1, nDim-1
            ! Last eight shape function derivative factors are identical to the first four, because the global position vector is defined using the shape functions of the top nodes (equation (7) in [4])
            DerivShapeFactors(j,i+nNodes/2) = DerivShapeFactors(j,i)
         END DO
      END DO
!
! ! Debugging
!       errout = erhandler('USERELEM',1008,2,
!      & 'ShapeFunc. Printing some values   
!      & %F, %F, %F, %F, %F.',
!      & (/eta, xi, nNodes*ONE, ShapeFactors(1),  
!      & DerivShapeFactors(1,1)/), ' ')
!
!     End of the subroutine
      RETURN
      END SUBROUTINE ShapeFunc
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE Jacobian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     Subroutine to calculate the Jacobian matrix of the element and to obtain the area scale factor for the current integration point
      SUBROUTINE Jacobian(nDim, nNodes, DerivShapeFactors, xRef, nlgeom,
     &  xCur, Theta, dA, NIntP, IntPWeights, iInt, jInt)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      INTEGER 
     &  nDim,                                      ! Number of dimensions of the problem (sub-in) 
     &  nNodes,                                    ! Number of nodes of the element (sub-in) 
     &  nlgeom,                                    ! Large deformation key [from nlgeom command] (sub-in)
     &  NIntP,                                     ! Number of integration points (sub-in)
     &  iInt,                                      ! Integration counter along natural coordinate eta (sub-in)
     &  jInt,                                      ! Integration counter along natural coordinate xi (sub-in)
     &  i,                                         ! Used in loops (inside)
     &  j,                                         ! Used in loops (inside)
     &  k                                          ! Used in loops (inside)
        DOUBLE PRECISION
     &  JacobianMatrix(nDim-1, nDim),              ! Jacobian matrix (derivative of the global coordinates wrt the local coordinates) (inside)
     &  DerivShapeFactors(nDim-1,nNodes),          ! Factors for the derivatives of the shape functions (sub-in)    
     &  xRef(nDim,nNodes),                         ! Nodal coordinates in initial configuration (sub-in)
     &  xCur(nDim,nNodes),                         ! Nodal coordinates in current configuration (sub-in)
     &  vsNorm,                                    ! Norm of the first tangential directional vector at the current integration point (inside)
     &  vnNorm,                                    ! Norm of the out-of-plane directional vector at the current integration point (inside)
     &  Theta(nDim, nDim),                         ! Matrix of directional vectors. Orthogonal rotation tensor relating the local and fixed coordinate systems (sub-out)
     &  dA,                                        ! Area scale factor for the current integration point (inside) (sub-out)
     &  IntPWeights(NIntP)                         ! Weights for the integration points (sub-in)
!    Next, the variables used in debugging are declared
        EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
        INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout                                     ! Output of erhandler (ANSYS own subroutine for displaying program errors) (inside)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!     Initialise the variables
      DO i = 1, nDim-1
        DO j = 1, nDim
            JacobianMatrix(i,j) = ZERO
        END DO
      END DO
!
!     Calculate the Jacobian matrix as the product of the shape functions by the global coordinates of the midsurface
      DO i = 1, nDim-1
        DO j = 1, nDim
            ! Loop over the number of nodes: J(i,j) = Sum(DerivShapeFactors(i, node)*xRef(j,node))
            DO k = 1, nNodes/2
                ! If there are no geometric non-linearities, employ the undeformed nodal coordinates
                IF (nlgeom == 0) THEN
                    JacobianMatrix(i,j) = JacobianMatrix(i,j) + 
     & DerivShapeFactors(i,k)*xRef(j,k)
                ! If there are geometric non-linearities, the midsurface coordinates will not coincide with the undeformed configuration
                ! They are thus calculated as the average between the bottom node and its top node counterpart current coordinates
                ELSE
                    JacobianMatrix(i,j) = JacobianMatrix(i,j) + 
     & DerivShapeFactors(i,k) * HALF*(xCur(j,k) + xCur(j,k+nNodes/2))
                END IF
            END DO
        END DO
      END DO
!
!     Obtain the orthogonal rotation tensor relating the local and fixed coordinate systems by calculating the directional vectors
      ! This is done as indicated in equations (9-11) in [1], such that Theta(1,:) is the first tangential directional vector vs, 
      ! Theta(2,:) is the second tangential directional vector vt, and Theta(3,:) is the normal directional vector vn
      ! It must be noted that eta and xi are defined the other way around in (7,8) in [1]
!     First, obtain the first tangential vector vs as the normalised dNi/deta. Equation (10) in [1] 
      vsNorm = SQRT(JacobianMatrix(1,1)**2 + JacobianMatrix(1,2)**2 + 
     & JacobianMatrix(1,3)**2)
      DO i = 1, 3
        Theta(1,i) = JacobianMatrix(1,i)/vsNorm
      END DO
!     Second, obtain the out-of-plane normal directional vector as the the cross product of veta and vxi, expressed through determinants (equation (9) in [1]) 
      Theta(3,1) = JacobianMatrix(1,2)*JacobianMatrix(2,3) - 
     & JacobianMatrix(2,2)*JacobianMatrix(1,3)
      Theta(3,2) = JacobianMatrix(2,1)*JacobianMatrix(1,3) - 
     & JacobianMatrix(1,1)*JacobianMatrix(2,3)
      Theta(3,3) = JacobianMatrix(1,1)*JacobianMatrix(2,2) - 
     & JacobianMatrix(2,1)*JacobianMatrix(1,2)
      ! Obtain the norm of the out-of-plane normal directional vector 
      vnNorm = SQRT(Theta(3,1)**2 + Theta(3,2)**2 + Theta(3,3)**2)
      ! Normalise the out-of-plane normal directional vector
      DO i = 1, 3
        Theta(3,i) = Theta(3,i)/vnNorm
      END DO
!     Third, obtain the second tangential vector as the cross product of the previous two. No need to normalise now. Equation (11) in [1]
      Theta(2,1) = Theta(3,2)*Theta(1,3) - Theta(3,3)*Theta(1,2)
      Theta(2,2) = Theta(3,3)*Theta(1,1) - Theta(3,1)*Theta(1,3)
      Theta(2,3) = Theta(3,1)*Theta(1,2) - Theta(3,2)*Theta(1,1)
!
!     Obtain the area scale factor for the current integration point using the norm of the out-of-plane vector and the weights of the current integration point
      dA = vnNorm*IntPWeights(iInt)*IntPWeights(jInt)
!
! !     Debugging
!       errout = erhandler('USERELEM',1009,2,
!      & 'Jacobian. Printing some values   
!      & %F, %F, %F, %F, %F, %F, %F.',
!      & (/JacobianMatrix(1,1), vsNorm, vnNorm,  
!      & Theta(3,1), Theta(3,2), Theta(3,3), dA/), ' ')
!
!     End of the subroutine
      RETURN
      END SUBROUTINE Jacobian
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE Displacements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     Subroutine to calculate the displacements in local coordinates, as well as the B matrix in global coordinates 
      SUBROUTINE Displacements(nDim, nUsrDof, nNodes, ShapeFactors,   
     &  Bglobal, Theta, DisplJump, NIntP, IntP, TotValDofs)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      INTEGER 
     &  nDim,                                      ! Number of dimensions of the problem (sub-in)
     &  nUsrDof,                                   ! Number of DOFs of this element (matrix and load vector size) (sub-in)
     &  nNodes,                                    ! Number of nodes of the element (sub-in)
     &  i,                                         ! Used in loops (inside)
     &  j,                                         ! Used in loops (inside)
     &  k,                                         ! Used in loops (inside)
     &  NIntP,                                     ! Number of integration points (sub-in)
     &  IntP                                       ! Integration point number (sub-in)
        DOUBLE PRECISION
     &  Blocal(nDim, nUsrDof),                     ! B-matrix in local coordinates (inside)
     &  ShapeFactors(nNodes),                      ! Shape function factors in local coordinates (sub-in)
     &  Bglobal(nDim, nUsrDof),                    ! B-matrix in global coordinates (sub-out)
     &  Theta(nDim, nDim),                         ! Matrix of directional vectors. Orthogonal rotation tensor relating the local and fixed coordinate systems (sub-in)
     &  DisplJump(NIntP**2, nDim),                 ! Displacement jump at integration point (sub-out)
     &  TotValDofs(nUsrDof)                        ! Total values of DOFs (displacements) from time 0 (sub-in)
!    Next, the variables used in debugging are declared
        EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
        INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout                                     ! Output of erhandler (ANSYS own subroutine for displaying program errors) (inside)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!     Initialise the variables
      DO i = 1, nDim
        DO j = 1, nUsrDof
            Blocal(i,j) = ZERO
            Bglobal(i,j) = ZERO
        END DO
        DisplJump(IntP, i) = ZERO
      END DO
!
!     Define the B-matrix in local coordinates using the shape function factors: B = [N1 0 0 N2 0 0 N3 0 0 ...; 0 N1 0 0 N2 0 0 N3 0 ...; 0 0 N1 0 0 N2 0 0 N3 ...]     
      DO i = 1, nDim
        k = 0
        DO j = i, nUsrDof, nDim
            k = k + 1
            Blocal(i,j) = ShapeFactors(k)
        END DO
      END DO
!     Define the B-matrix in global coordinates as the product of the orthogonal rotation tensor relating the local and fixed coordinate systems Theta and the previous B-matrix     
      DO i = 1, nDim
        DO j = 1, nUsrDof
            DO k = 1, nDim
                Bglobal(i,j) = Bglobal(i,j) + Theta(i, k)*Blocal(k,j)
            END DO
        END DO
      END DO
!
!     Obtain the displacement jump in local coordinates
      ! This is equation (66) in [1], where the calculation of the global B-matrix is also indicated
      DO i = 1, nDim
        DO j = 1, nUsrDof
            DisplJump(IntP, i) = DisplJump(IntP, i) + 
     & Bglobal(i,j)*TotValDofs(j)
        END DO
      END DO
!
! !     Debugging
!       errout = erhandler('USERELEM',1010,2,
!      & 'Displacements. Printing some values   
!      & %F, %F, %F. 
!      & Nodal displacements: %F, %F, %F, %F, %F
!      & %F, %F, %F, %F, %F, %F, %F, %F, %F, %F
!      & %F, %F, %F, %F, %F, %F, %F, %F, %F.',
!      & (/Blocal(1,1), Bglobal(1,1), DisplJump(1, 1),
!      & TotValDofs(1), TotValDofs(2), TotValDofs(3),
!      & TotValDofs(4), TotValDofs(5), TotValDofs(6),
!      & TotValDofs(7), TotValDofs(8), TotValDofs(9),
!      & TotValDofs(10), TotValDofs(11), TotValDofs(12),
!      & TotValDofs(13), TotValDofs(14), TotValDofs(15),
!      & TotValDofs(16), TotValDofs(17), TotValDofs(18),
!      & TotValDofs(19), TotValDofs(20), TotValDofs(21),
!      & TotValDofs(22), TotValDofs(23), TotValDofs(24)/), ' ')
!
!     End of the subroutine
      RETURN
      END SUBROUTINE Displacements
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE Damage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     Subroutine to calculate the stiffness damage variable and thus the stiffness constitutive matrix and parameter KH for the calculation of the tangent matrix
      ! The subroutine first characterises the cohesive law based on the current displacement-based mode mixity
      ! Then, the subroutine Rates is called to obtain the derivative of the energy damage variable with respect to the crack growth 
      ! Next, the subroutines DamageStatic or DamageFatigue are called to calculate the evolution of the stiffness damage variable Dk under static or fatigue loading
      ! For fatigue initiation, this subroutine also calls InitiationVariables to calculate the initiation damage variable and other fatigue initiation variables
      ! Finally, the stiffness constitutive matrix and parameter KH for the calculation of the tangent matrix are obtained
      SUBROUTINE Damage(nDim, NIntP, IntP, DisplJump, tauI, tauII,  
     &  KStiffness, BKeta, GIc, GIIc, nSaveVars, saveVars, ldstep,    
     &  Loading, DkNew, timval, GI, GII, PropCrack, cParisI, 
     &  cParisII, cParism, mParisI, mParisII, mParism, DeltaN, 
     &  kfstps, ddelta3da, ddelta1da, sm, SeverityInit, DiNew, sI, 
     &  sII, DeltaNi, InitiatedCrack, ConstStiff, KH, keyHisUpd,   
     &  ieqitr, isubst, elId)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      INTEGER
     &  nDim,                                      ! Number of dimensions of the problem (sub-in)
     &  NIntP,                                     ! Number of integration points (sub-in)
     &  IntP,                                      ! Integration point number (sub-in)
     &  nSaveVars,                                 ! Number of saved variables (sub-in)
     &  ldstep,                                    ! Current load step number (sub-in)
     &  Loading(NIntP**2),                         ! Vector of flags indicating the loading or not of the cohesive element (sub-out)
     &  PropCrack,                                 ! Variable indicating whether a crack should go through fatigue propagation damage (1) or fatigue initiation damage (0)         
     &  kfstps,                                    ! Key for the first iteration of first substep of the first load step (sub-in)
     &  InitiatedCrack,                            ! Variable indicating whether this element belongs to an initiated crack (sub-in)
     &  keyHisUpd,                                 ! Key to indicate if history-dependent variables need to be updated, since the iteration is already converged (sub-in)     
     &  i,                                         ! Used in loops (inside)
     &  j                                          ! Used in loops (inside)     
        DOUBLE PRECISION
     &  DisplJump(NIntP**2, nDim),                 ! Displacement jump at integration point (sub-in)
     &  deltaI,                                    ! Mode I displacement jump (inside)
     &  deltaS,                                    ! Shear (modes II and III) displacement jump (inside)
     &  delta,                                     ! Equivalent one-dimensional opening displacement (inside)
     &  beta,                                      ! Local opening displacement-based mode mixity (inside)
     &  B,                                         ! Energy-based mode mixity (inside)
     &  deltaI0,                                   ! Mode I onset separation (inside)
     &  deltaS0,                                   ! Shear (modes II and III) onset separation (inside)
     &  tauI,                                      ! Mode I interface strength (sub-in)
     &  tauII,                                     ! Mode II interface strength (sub-in)
     &  KStiffness,                                ! Penalty stiffness. Initial stiffness of undamaged cohesive element (sub-in)
     &  delta0,                                    ! Mixed mode onset equivalent displacement (inside)
     &  BKeta,                                     ! Benzeggagh–Kenane criterion coefficient (sub-in)
     &  Gc,                                        ! Mixed mode critical energy release rate (inside)
     &  GIc,                                       ! Mode I critical energy release rate / fracture toughness (sub-in)
     &  GIIc,                                      ! Mode II critical energy release rate / fracture toughness (sub-in)
     &  tau0,                                      ! Mixed mode equivalent interface strength (inside)
     &  deltac,                                    ! Mixed mode final equivalent displacement (inside)
     &  saveVars(nSaveVars),                       ! User saved variables (sub-inout)
     &  DkOld,                                     ! Stiffness damage variable from the previous substep (inside)
     &  deltath,                                   ! Threshold damage variable (inside)
     &  DkNew,                                     ! Stiffness damage variable in the current substep (sub-out)
     &  timval,                                    ! Current time value (sub-in)
     &  TimeOld,                                   ! Time in the previous substep (inside)
     &  delta0Old,                                 ! Mixed mode onset equivalent displacement in the previous increment (inside)
     &  deltacOld,                                 ! Mixed mode final equivalent displacement in the previous increment (inside)
     &  wtotOld,                                   ! Total dissipated specific energy in the previous increment (inside)
     &  wtot,                                      ! Total dissipated specific energy in the current increment (inside)
     &  GcOld,                                     ! Mixed mode critical energy release rate in the previous increment (inside)
     &  GI,                                        ! Mode I energy release rate obtained from the J-integral calculation (sub-in)
     &  GII,                                       ! Mode II energy release rate obtained from the J-integral calculation (sub-in)
     &  cParisI,                                   ! Parameter c in Paris law for mode I (sub-in)
     &  cParisII,                                  ! Parameter c in Paris law for mode II (sub-in)
     &  cParism,                                   ! Parameter c in Paris law for fitting of mode interaction (sub-in)
     &  mParisI,                                   ! Parameter m in Paris law for mode I (sub-in)
     &  mParisII,                                  ! Parameter m in Paris law for mode II (sub-in)
     &  mParism,                                   ! Parameter m in Paris law for fitting of mode interaction (sub-in)     
     &  ddelta3da,                                 ! Derivative of the opening displacement with respect to the crack growth a (sub-in)
     &  ddelta1da,                                 ! Derivative of the mode II displacement with respect to the crack growth a (sub-in)
     &  DeltaN,                                    ! Increment of cycles used in the calculations (sub-in)
     &  delta3,                                    ! Opening displacement jump at the integration point (inside)   
     &  dDedaOld,                                  ! Derivative of the energy damage variable with respect to the crack growth from the previous substep (inside)   
     &  dDeda,                                     ! Derivative of the energy damage variable with respect to the crack growth (inside)
     &  Nloc,                                      ! Number of cycles in true local time (inside)
     &  Ncor,                                      ! Number of cycles in corrected local time (inside)
     &  Ns,                                        ! Number of cycles related to the unintended evolution of quasi-static damage (inside)
     &  sm,                                        ! Mode mixity S-N curve slope (sub-out)
     &  SeverityInit,                              ! Severity: ratio of current traction over static strength (sub-out)
     &  sI,                                        ! Mode I initiation S-N curve slope (sub-in)
     &  sII,                                       ! Mode II initiation S-N curve slope (sub-in)
     &  DeltaNi,                                   ! Increment of cycles based on initiation for next substep (sub-out)
     &  DiNew,                                     ! Initiation damage variable in the current cycle increment (sub-out)
     &  DiOld,                                     ! Initiation damage variable in the previous cycle increment (sub-in)
     &  ConstStiff(NIntP**2, nDim, nDim),          ! Constitutive stiffness matrix (sub-out)
     &  KH(NIntP**2)                               ! Parameter used in the calculation of the constitutive tangent stiffness matrix (sub-out)
!    Next, the variables used in debugging are declared. Including some variables used to identify the instant when the message is printed out
        EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
        INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout,                                    ! Output of erhandler (ANSYS own subroutine for displaying program errors) (inside)   
     &  ieqitr,                                    ! Current equilibrium iteration number (sub-in)                                    
     &  isubst,                                    ! Current substep number (sub-in)
     &  elId                                       ! Element number (sub-in)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!     Initialise the variables
      DO i = 1, nDim
        DO j = 1, nDim
            ConstStiff(IntP, i, j) = ZERO
        END DO
      END DO
      Loading(IntP) = 0
!
!     Read variables from the previous substep for this element
      DkOld = saveVars(IntP)
      delta0Old = saveVars((NIntP**2)*1 + IntP)
      deltacOld = saveVars((NIntP**2)*2 + IntP)
      wtotOld = saveVars((NIntP**2)*3 + IntP)
      TimeOld = saveVars((NIntP**2)*4 + IntP) 
      GcOld = saveVars((NIntP**2)*5 + IntP) 
      dDedaOld = saveVars((NIntP**2)*6 + IntP)
      Nloc = saveVars((NIntP**2)*7 + IntP)
      Ncor = saveVars((NIntP**2)*8 + IntP)
      Ns = saveVars((NIntP**2)*9 + IntP)
      DiOld = saveVars((NIntP**2)*11 + IntP)
!
!     Characterise the cohesive law
!     Calculate the equivalent opening displacement. 
      ! Equations from (8) in [2] are employed
      deltaI = HALF*(DisplJump(IntP, 3) + ABS(DisplJump(IntP, 3)))
      deltaS = SQRT(DisplJump(IntP, 1)**2 + DisplJump(IntP, 2)**2)
      delta = SQRT(deltaI**2 + deltaS**2)
!     Calculate mode mixity and mixed mode onset and final displacements with B-KStiffness criterion. 
      ! Equation (9) from [2]
      beta = deltaS / (deltaI + deltaS)
      ! In order to avoid beta becoming NAN if the displacement jumps are 0 (and thus both deltaS and deltaI), impose mode II for very small deltaI, no matter deltaII
      IF (deltaI < 1.0D-16) THEN
            beta = ONE 
      END IF
      ! Equation (11) from [2]
      B = beta**2 / (TWO*beta**2 - TWO*beta + ONE)
      ! Obtain the mode I and shear mode onset displacements. Equation (10) from [2]  
      deltaI0 = tauI/KStiffness
      deltaS0 = tauII/KStiffness
      ! Obtain the the mixed mode onset displacement. Equation (45) from [1]  (with different nomenclature)
      delta0 = SQRT(deltaI0**2 + (deltaS0**2 - deltaI0**2)*B**BKeta)
      ! Obtain the the mixed mode critical energy release rate. Equation (11) from [2] 
      Gc = GIc + (GIIc - GIc)*B**BKeta
      ! Obtain the mixed mode onset strength and final displacement. Equations (10) from [2]
      tau0 = delta0*KStiffness
      deltac = TWO*Gc/tau0
!     Calculate the damage threshold: damage state based on the stiffness damage variable from the previous substep 
      ! Equation (14) from [2] for the opening displacement associated to a damage status
      deltath =  delta0*deltac / (deltac - DkOld*(deltac - delta0))
      ! If the stiffness damage variable from the previous substep is equal or greater than 1, then the threshold value value must be deltac
      IF (DkOld .GE. ONE) THEN
        deltath = deltac
      END IF
      ! Generate an error message if the threshold value is greater than the final displacement
      IF (deltath > deltac) THEN
            errout = erhandler('USERELEM',5003,3,
     & 'ERROR. deltath %F is greater than deltac %F.
     & Other values are delta0 %F, DkOld %F,
     & deltaI0 %F, deltaS0 %F, Gc %F, tau0 %F',
     & (/deltath, deltac, delta0, DkOld,
     & deltaI0, deltaS0, Gc, tau0/), ' ')
      END IF
!
!     No matter if the loading is quasi-static or fatigue, the derivative of the energy damage variable with respect to the crack growth must be calculated
      ! In static this is necessary so that this derivative can be used in the first substep of the following fatigue-loading substep 
      ! The subroutine Rates is employed for this calculation
      CALL Rates(KStiffness, deltac, delta, delta0, Gc, tau0, 
     & BKeta, B, tauI, tauII, GIc, GIIc, beta, DisplJump(IntP, 3), 
     & ddelta3da, ddelta1da, deltaI, deltaS, dDeda, ieqitr, isubst,
     & elId)
!
!     Calculate the stiffness damage variable Dk depending on the type of loading
!     In the first load step, quasi-static loading is applied. Also in the case that the element does not belong to a fatigue 'propagation crack'
      IF (ldstep == 1 .OR. PropCrack == 0) THEN
            ! Call subroutine DamageStatic
            CALL DamageStatic(delta, deltath, NIntP, IntP, Loading, 
     & ieqitr, DkNew, deltac, delta0, DkOld, kfstps, KStiffness, 
     & delta0Old, deltacOld, tau0, wtotOld, wtot, isubst, elId)
!     In the subsequent load steps and if the element is part of a 'propagation crack', fatigue loading is applied
      ELSEIF (ldstep > 1 .AND. PropCrack == 1) THEN 
            ! Call subroutine DamageFatigue
            CALL DamageFatigue(DkOld, delta0Old, deltacOld, KStiffness,
     & GcOld, GI, GII, cParisI, cParisII, cParism, mParisI, mParisII, 
     & mParism, GIc, GIIc, BKeta, DeltaN, delta, delta0, dDeda,  
     & dDedaOld, Nloc, Ncor, Ns, deltac, Gc, DkNew, tau0, wtotOld, wtot,  
     & ieqitr, isubst, elId)
      END IF
!
!     If the element belongs to an 'initiation crack', it is fully damaged, i.e., the stiffness damage variable Dk is 1
      IF (InitiatedCrack == 1) THEN
            DkNew = ONE
            Loading = 0
      END IF
!
!     Calculate fatigue initiation variables through subroutine InitiationVariables
      ! Done after after the first static substep, since in the first substep some variables are being zeroed
      IF ((ldstep == 1 .AND. isubst > 1) .OR. ldstep > 1) THEN
            CALL InitiationVariables(sI, sII, mParisI, mParisII, 
     & deltaI, delta, deltaS, KStiffness, DkNew, sm, SeverityInit,  
     & tau0, DiNew, DiOld, DeltaN, DeltaNi, 
     & ieqitr, isubst, elId)
      END IF 
!
!     If the iteration has converged, then update the save variables to be accessed only the first time the next substep is run for this element
      IF (keyHisUpd == 1) THEN
            saveVars(IntP) = DkNew
            saveVars((NIntP**2)*1 + IntP) = delta0
            saveVars((NIntP**2)*2 + IntP) = deltac
            saveVars((NIntP**2)*3 + IntP) = wtot
            saveVars((NIntP**2)*4 + IntP) = timval
            saveVars((NIntP**2)*5 + IntP) = Gc
            saveVars((NIntP**2)*6 + IntP) = dDeda
            saveVars((NIntP**2)*7 + IntP) = Nloc
            saveVars((NIntP**2)*8 + IntP) = Ncor
            saveVars((NIntP**2)*9 + IntP) = Ns
            saveVars((NIntP**2)*11 + IntP) = DiNew
      END IF
!
!     Calculate the constitutive stiffness matrix
      ! Equation (47) from [1]. Everything is zero, except for the diagonal, which is (1-DkNew)*KStiffness
      DO i = 1, nDim
        ConstStiff(IntP, i, i) = (ONE-DkNew)*KStiffness 
      END DO
      ! With the exception of penetration, if the normal displacement is negative, then it is KStiffness to avoid such penetration
      IF (DisplJump(IntP, 3) < ZERO) THEN 
        ConstStiff(IntP, 3, 3) = KStiffness
      END IF
!
!     Calculate the KH parameter used in the calculation of the constitutive tangent stiffness if the cohesive element is being loaded
      ! Equation (53) from [1], multiplied by K, as done in equation (52)
      IF (Loading(IntP) == 1) THEN
        KH(IntP) = KStiffness*deltac*delta0 / 
     & (deltath**3*(deltac - delta0))
      END IF
!
! !     Debugging
!       errout = erhandler('USERELEM',1011,2,
!      & 'Damage. Printing some values:
!      & Substep %F Iteration %F Element %F Time %F
!      & KeyConverged %F.
!      & Deltas: %F, %F, %F, %F, %F, %F.
!      & Origin of delta in pure mode I: %F, %F.
!      & Times: %F, %F, %F, %F.
!      & Const Stiff Matrix: %F, %F, %F, %F, %F, %F, %F.
!      & KH: %F, %F.',
!      & (/ isubst*1.D0, ieqitr*1.D0, elId*1.D0, timval*1.D0, 
!      & keyHisUpd*1.D0,
!      & deltac, delta0, delta, deltath, DkOld, DkNew,
!      & deltaI, DisplJump(IntP, 3), 
!      & timval, TimeOld,    
!      & saveVars((NIntP**2)*4 + IntP), saveVars((NIntP**2)*5 + IntP),
!      & ConstStiff(1, 1, 1), ConstStiff(1, 1, 1)/1e6,
!      & ConstStiff(1, 2, 2)/1e6, ConstStiff(1, 3, 3)/1e6, 
!      & (ONE-DkNew)*KStiffness, DkNew, KStiffness, 
!      & KH(IntP)/1e6, Loading(IntP)*ONE/), 
!      & ' ')
!
!     End of the subroutine
      RETURN
      END SUBROUTINE Damage
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE Rates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Subroutine to calculate the derivative of the energy damage variable with respect to the crack growth (dDeda)
      SUBROUTINE Rates(KStiffness, deltac, delta, delta0, Gc, tau0, 
     & BKeta, B, tauI, tauII, GIc, GIIc, beta, delta3, ddelta3da, 
     & ddelta1da, deltaI, deltaS, dDeda, ieqitr, isubst, elId)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      DOUBLE PRECISION
     &  wr,                                        ! Specific remaining ability to do non-conservative work (inside)
     &  KStiffness,                                ! Penalty stiffness. Initial stiffness of undamaged cohesive element (sub-in)
     &  deltac,                                    ! Mixed mode final equivalent displacement (sub-in)
     &  delta,                                     ! Equivalent one-dimensional opening displacement (sub-in)
     &  delta0,                                    ! Mixed mode onset equivalent displacement (sub-in)
     &  dDedwr,                                    ! Derivative of the energy damage variable wrt the specific remaining ability to do non-conservative work (inside)
     &  Gc,                                        ! Mixed mode local critical energy release rate (sub-in)
     &  dwrddelta,                                 ! Derivative of the specific remaining ability to do non-conservative work wrt the equivalent one-dimensional displacement (inside)
     &  dwrddelta0,                                ! Derivative of the specific remaining ability to do non-conservative work wrt the onset equivalent displacement (inside)
     &  ddelta0dtau0,                              ! Derivative of the onset equivalent displacement wrt the onset equivalent traction (inside)
     &  dwrddeltac,                                ! Derivative of the specific remaining ability to do non-conservative work wrt the final equivalent displacement (inside)
     &  ddeltacdtau0,                              ! Derivative of the final equivalent displacement wrt the onset equivalent traction (inside)   
     &  tau0,                                      ! Mixed mode equivalent interface strength (sub-in)       
     &  dtau0dB,                                   ! Derivative of the onset equivalent traction wrt the energy-based mode mixity (inside)
     &  BKeta,                                     ! Benzeggagh–Kenane criterion coefficient (sub-in)
     &  B,                                         ! Local energy-based mode mixity (sub-in)
     &  tauI,                                      ! Mode I interface strength (sub-in)
     &  tauII,                                     ! Mode II interface strength (sub-in)
     &  ddeltacdGc,                                ! Derivative of the final equivalent displacement wrt the local critical energy release rate (inside)
     &  dDedGc,                                    ! Derivative of the energy damage varaible wrt the local critical energy release rate (inside)
     &  dGcdB,                                     ! Derivative of the local critical energy release rate wrt local energy-based mode mixity (inside)
     &  GIc,                                       ! Mode I critical energy release rate / fracture toughness (sub-in)
     &  GIIc,                                      ! Mode II critical energy release rate / fracture toughness (sub-in)
     &  dBdbeta,                                   ! Derivative of the local energy-based mode mixity wrt local opening displacement-based mode mixity (inside)
     &  beta,                                      ! Local opening displacement-based mode mixity (sub-in)
     &  delta3,                                    ! Opening displacement jump at the integration point (sub-in)
     &  ddeltada,                                  ! Derivative of the equivalent one-dimensional displacement wrt crack growth (inside)
     &  dbetada,                                   ! Derivative of the local opening displacement-based mode mixity wrt crack growth (inside)
     &  ddelta3da,                                 ! Derivative of the opening displacement with respect to the crack growth a (sub-in)
     &  ddelta1da,                                 ! Derivative of the mode II displacement with respect to the crack growth a (sub-in)
     &  deltaI,                                    ! Mode I displacement jump (sub-in)
     &  deltaS,                                    ! Shear (modes II and III) displacement jump (sub-in)
     &  dDeda                                      ! Derivative of the energy damage variable with respect to the crack growth (sub-out)
!    Next, the variables used in debugging are declared
        EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
        INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout,                                    ! Output of erhandler (ANSYS own subroutine for displaying program errors) (inside)   
     &  ieqitr,                                    ! Current equilibrium iteration number (sub-in)                                    
     &  isubst,                                    ! Current substep number (sub-in)
     &  elId                                       ! Element number (sub-in)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!     Firstly, calculate the specific remaining ability to do non-conservative work wr. Equation (21) from [1]
      wr = KStiffness*delta0*deltac*(deltac-delta) /
     & (TWO*(deltac-delta0))
!     Secondly, calculate the partial derivatives from equations (24) from [1]
      dDedwr = -ONE/Gc
      dwrddelta = KStiffness*delta0*deltac /
     & (TWO*(delta0-deltac))
      dwrddelta0 = KStiffness*(deltac-delta)*deltac**2 /
     & (TWO*(deltac-delta0)**2)
      ddelta0dtau0 = ONE/KStiffness
      dwrddeltac = KStiffness*delta0*(delta*delta0 
     & + deltac**2 - TWO*deltac*delta0) /
     & (TWO*(deltac-delta0)**2)
      ddeltacdtau0 = -TWO*Gc / tau0**2
      dtau0dB = BKeta*(tauII**2-tauI**2)*B**(BKeta-ONE) / 
     & (TWO*tau0)
      ddeltacdGc = TWO/tau0
      dDedGc = wr/Gc**2
      dGcdB = BKeta*(GIIc-GIc)*B**(BKeta-ONE)  
      dBdbeta = (-TWO*beta**2 + TWO*beta)
     & / (TWO*beta**2 - TWO*beta + ONE)**2
!     Thirdly, calculate the derivatives of the equivalent displacement and displacement-based mode mixity wrt the crack growth
      ! Equation (36) from [2] for ddeltada and equation (37) for dbetada
      ! In equation (37), ddelta1da is in reality ddeltaSdx, so the absolute value of the slope ddelta1da must be used
      IF (delta3 > ZERO) THEN
            ddeltada = deltaI/delta*ddelta3da +
     & deltaS/delta*ABS(ddelta1da)
            dbetada = deltaI/((deltaS+deltaI)**2)*
     & abs(ddelta1da) -
     & deltaS/((deltaS+deltaI)**2)*ddelta3da
      ELSE
            ddeltada = deltaS/delta*ABS(ddelta1da)
            dbetada = ZERO
      END IF
      ! Avoid undesired values if the equivalent displacement is 0
      IF (delta < 1.0D-10) THEN
            ddeltada = ZERO
            dbetada = ZERO
      END IF
!     Finally, obtain the derivative of the energy damage variable with respect to the crack growth
      ! Equation (22) from [2]
      dDeda = ( ( (dwrddelta0*ddelta0dtau0 + 
     & dwrddeltac*ddeltacdtau0)*dDedwr*dtau0dB + 
     & (dDedwr*dwrddeltac*ddeltacdGc + dDedGc)*dGcdB ) 
     & * dBdbeta * dbetada 
     & + (dDedwr*dwrddelta) * ddeltada)
!
! ! Debugging.
!       errout = erhandler('USERELEM',1012,2,
!      & 'Warning in dDeda %F
!      & Substep %F, iteration %F and element %F
!      & Derivatives: %F, %F, %F, %F.
!      & In detail: %F, %F, %F, %F, %F.
!      & Delta3: %F, %F, %F, %F, %F. ',
!      & (/dDeda,
!      & isubst*ONE, ieqitr*ONE, elId*ONE,
!      & dbetada, ddeltada,
!      & ( (dwrddelta0*ddelta0dtau0 + 
!      & dwrddeltac*ddeltacdtau0)*dDedwr*dtau0dB + 
!      & (dDedwr*dwrddeltac*ddeltacdGc + dDedGc)*dGcdB ) 
!      & * dBdbeta * dbetada,
!      & (dDedwr*dwrddelta) * ddeltada,
!      & dDedwr, dwrddelta, Gc, delta0, deltac,
!      & delta3, deltaI, delta, ddelta3da, 
!      & ddelta1da, deltaS/), ' ')
!
      RETURN
      END SUBROUTINE Rates
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE DamageStatic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Subroutine to calculate damage under static loading. 
      ! It also obtains the total dissipated specific energy, which is necessary for fatigue propagation
      SUBROUTINE DamageStatic(delta, deltath, NIntP, IntP, Loading, 
     & ieqitr, DkNew, deltac, delta0, DkOld, kfstps, KStiffness, 
     & delta0Old, deltacOld, tau0, wtotOld, wtot, isubst, elId)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      INTEGER
     &  NIntP,                                     ! Number of integration points (sub-in)
     &  IntP,                                      ! Integration point number (sub-in)
     &  Loading(NIntP**2),                         ! Vector of flags indicating the loading or not of the cohesive element (sub-out) 
     &  ieqitr,                                    ! Current equilibrium iteration number (sub-in)               
     &  kfstps                                     ! Key for the first iteration of first substep of the first load step (sub-in)      
        DOUBLE PRECISION
     &  delta,                                     ! Equivalent one-dimensional opening displacement (sub-in)
     &  deltath,                                   ! Threshold damage variable (sub-inout)
     &  DkNew,                                     ! Stiffness damage variable in the current substep (sub-out)
     &  deltac,                                    ! Mixed mode final equivalent displacement (sub-in)
     &  delta0,                                    ! Mixed mode onset equivalent displacement (sub-in)
     &  DkOld,                                     ! Stiffness damage variable from the previous substep (sub-in)
     &  deltathNew,                                ! Threshold damage variable in the current increment (inside)
     &  tauNew,                                    ! Mixed mode equivalent interface strength in the current increment (inside)
     &  KStiffness,                                ! Penalty stiffness. Initial stiffness of undamaged cohesive element (sub-in)
     &  deltathOld,                                ! Threshold damage variable in the previous increment (inside)
     &  delta0Old,                                 ! Mixed mode onset equivalent displacement in the previous increment (sub-in)
     &  deltacOld,                                 ! Mixed mode final equivalent displacement in the previous increment (sub-in)
     &  tauOld,                                    ! Mixed mode equivalent interface strength in the previous increment (inside)
     &  tau0,                                      ! Mixed mode equivalent interface strength (sub-in)
     &  wtotOld,                                   ! Total dissipated specific energy in the previous increment (sub-in)
     &  wtot                                       ! Total dissipated specific energy in the current increment (sub-out)   
!    Next, the variables used in debugging are declared. Including some variables used to identify the instant when the message is printed out
        EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
        INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout,                                    ! Output of erhandler (ANSYS own subroutine for displaying program errors) (inside)  
     &  isubst,                                    ! Current substep number (sub-in)
     &  elId                                       ! Element number (sub-in)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!     Calculate the new stiffness damage variable
      ! If delta is greater than deltath, the cohesive element is being loaded
      IF (delta > deltath) THEN ! In some previous versions of the subroutine, the condition was IF (delta > deltath .AND. delta < deltac) THEN
        ! The threshold is then the new equivalent opening displacement
        deltath = delta
        Loading(IntP) = 1
        ! If convergence was not achieved in the first iteration and full damage is expected, then use the secant stiffness matrix
        IF (ieqitr > 1 .AND. delta > deltac) THEN
            Loading(IntP) = 0
        END IF
      ! If not, the element is not being loaded and the threshold is not modified
      ELSE
        Loading(IntP) = 0
      END IF
      ! Calculate the new stiffness damage variable. Equation (13) from [2]
      DkNew = deltac*(deltath - delta0) / (deltath*(deltac - delta0))
      ! If DkNew is negative, then it must be 0. And if it is greater than 1, then it must be 1. Equation (12) from [2]
      IF (DkNew < ZERO) THEN
        DkNew = ZERO
      ELSEIF (DkNew > ONE) THEN
        DkNew = ONE
      END IF
!
!     Calculate the total dissipated specific energy. The total dissipated energy is the area in figure 4 of [2]. This is only needed for the fatigue damage section 
      ! If no damage has occured
      IF (DkNew < 0.1d-15) THEN
            ! In this case, from the area, the energy is a triangle of base delta and height KStiffness*delta
            wtot = HALF*KStiffness*delta**2
      ! If damage has occured in this substep
      ELSE
            ! Obtain the current damage state. Equation (14) from [2]
            deltathNew =  delta0*deltac / 
     & (deltac - DkNew*(deltac - delta0))
            ! Calculate the current traction. Equation (8) from [2]
            tauNew = (ONE - DkNew)*KStiffness*delta
            ! If it is not the first substep and damage has occured in the previous one, calculate the displacement and traction in the previous substep
            IF (kfstps .NE. 1 .AND. DkOld > ZERO) THEN
                  deltathOld =  delta0Old*deltacOld / 
     & (deltacOld - DkOld*(deltacOld - delta0Old))
                  tauOld = (ONE - DkOld)*KStiffness*deltathOld
            ! If it is the first susbtep or damage has not occured in the previous substep, the displacement and traction in the previous substep are the onset ones
            ELSE
                  deltathOld = delta0
                  tauOld = tau0
                  wtotOld = HALF*tau0*delta0
            END IF
            ! The dissipated energy is calculated in this case as the sum of the previously dissipated and the newly dissipated
            ! The newly dissipated specific energy is obtained by solving equation (17) from [2] via the trapezoidal rule
            ! With respect to figure 4, for the area of the trapezoid, the height is deltathNew-deltathOld and the bases are tauOld and tauNew   
            wtot = wtotOld + HALF*(tauOld + tauNew) * 
     & MAX(ZERO, deltathNew - deltathOld)
      END IF
      ! Zero if it is negative
      IF (wtot < ZERO) THEN
            wtot = ZERO
      END IF
!
! !     Debugging
!       errout = erhandler('USERELEM',1013,2,
!      & 'Static-loading damage. Printing some values:
!      & Substep %F Element %F
!      & Deltas: %F, %F, %F, %F, %F, %F.
!      & Calculation Dk: %F',
!      & (/ isubst*1.D0, elId*1.D0, 
!      & deltac, delta0, delta, deltath, DkOld, DkNew,
!      & deltac*(deltath - delta0) / (deltath*(deltac - delta0))/), 
!      & ' ')
!
      RETURN
      END SUBROUTINE DamageStatic
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE DamageFatigue %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Subroutine to calculate damage under fatigue loading
      ! First, the previous energy damage variable is calculated from the previous stiffness damage variable
      ! Second, the crack growth rate da/dN is calculated using the Paris' law and the J-integral calculations
      ! Third, the new stiffness damage variable Dk is calculated. A åpredictor-corrector scheme is used to account for unintended quasi-static damage
      ! Last, the dissipated specific energy is calculated to apply the energy full-damage criterion
      SUBROUTINE DamageFatigue(DkOld, delta0Old, deltacOld, KStiffness,
     & GcOld, GI, GII, cParisI, cParisII, cParism, mParisI, mParisII, 
     & mParism, GIc, GIIc, BKeta, DeltaN, delta, delta0, dDeda,  
     & dDedaOld, Nloc, Ncor, Ns, deltac, Gc, DkNew, tau0, wtotOld, wtot,  
     & ieqitr, isubst, elId)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine     
        DOUBLE PRECISION
     &  DkOld,                                     ! Stiffness damage variable from the previous substep (sub-in)
     &  DeOld,                                     ! Energy damage variable based on the stiffness damage variable from the previous substep (inside)
     &  deltaDOld,                                 ! One-dimensional opening displacement associated to the damage state from the previous substep (inside)
     &  delta0Old,                                 ! Mixed mode onset equivalent displacement in the previous increment (sub-in)
     &  deltacOld,                                 ! Mixed mode final equivalent displacement in the previous increment (sub-in)
     &  KStiffness,                                ! Penalty stiffness. Initial stiffness of undamaged cohesive element (sub-in)
     &  GcOld,                                     ! Mixed mode critical energy release rate in the previous increment (sub-in)
     &  G,                                         ! Total energy release rate obtained from the J-integral calculation (inside)
     &  GI,                                        ! Mode I energy release rate obtained from the J-integral calculation (sub-in)
     &  GII,                                       ! Mode II energy release rate obtained from the J-integral calculation (sub-in)
     &  modemixity,                                ! Mode mixity based on the J-integral calculation (inside)
     &  cParis,                                    ! Parameter c in Paris law for a certain mode mixity (inside)
     &  cParisI,                                   ! Parameter c in Paris law for mode I (sub-in)
     &  cParisII,                                  ! Parameter c in Paris law for mode II (sub-in)
     &  cParism,                                   ! Parameter c in Paris law for fitting of mode interaction (sub-in)
     &  mParis,                                    ! Parameter m in Paris law for a certain mode mixity (inside)
     &  mParisI,                                   ! Parameter m in Paris law for mode I (sub-in)
     &  mParisII,                                  ! Parameter m in Paris law for mode II (sub-in)
     &  mParism,                                   ! Parameter m in Paris law for fitting of mode interaction (sub-in)
     &  GcGlobal,                                  ! Mixed mode global critical energy release rate (inside)    
     &  GIc,                                       ! Mode I critical energy release rate / fracture toughness (sub-in)
     &  GIIc,                                      ! Mode II critical energy release rate / fracture toughness (sub-in)
     &  BKeta,                                     ! Benzeggagh–Kenane criterion coefficient (sub-in)
     &  dadN,                                      ! Crack growth rate: derivative of the crack growth wrt the number of cycles (inside)
     &  DeltaN,                                    ! Increment of cycles used in the calculations (sub-in)
     &  delta,                                     ! Equivalent one-dimensional opening displacement (sub-in)
     &  delta0,                                    ! Mixed mode onset equivalent displacement (sub-in)   
     &  dDeda,                                     ! Derivative of the energy damage variable with respect to the crack growth (sub-in)
     &  dDedNNew,                                  ! Derivative of the energy damage variable with respect to the number of cycles in the current substep (inside)    
     &  dDedaOld,                                  ! Derivative of the energy damage variable with respect to the crack growth from the previous substep (sub-in)
     &  dDedNOld,                                  ! Derivative of the energy damage variable with respect to the number of cycles in the previous substep (inside)
     &  DeltaNcor,                                 ! Increment of cycles in corrected local time (inside)
     &  Nloc,                                      ! Number of cycles in true local time (sub-inout)
     &  Ncor,                                      ! Number of cycles in corrected local time (sub-inout)
     &  Ns,                                        ! Number of cycles related to the unintended evolution of quasi-static damage (sub-inout)
     &  DeNew,                                     ! Energy damage variable in the current substep (inside)
     &  deltath,                                   ! Threshold damage variable (inside)
     &  deltac,                                    ! Mixed mode final equivalent displacement (sub-in)
     &  Gc,                                        ! Mixed mode local critical energy release rate (sub-in)  
     &  DkNew,                                     ! Stiffness damage variable in the current substep (sub-out)
     &  DefatigueNew,                              ! Energy damage variable due solely to fatigue when there is unintended quasi-static damage (inside)
     &  DeltaDes,                                  ! Unintended quasi-static damage (inside)
     &  tauNew,                                    ! Mixed mode equivalent interface strength in the current increment (inside)
     &  deltathNew,                                ! Threshold damage variable in the current increment (inside)
     &  deltathOld,                                ! Threshold damage variable in the previous increment (inside)
     &  tauOld,                                    ! Mixed mode equivalent interface strength in the previous increment (inside)
     &  tau0,                                      ! Mixed mode equivalent interface strength (sub-in)
     &  wtotOld,                                   ! Total dissipated specific energy in the previous increment (sub-in)
     &  wtot                                       ! Total dissipated specific energy in the current increment (sub-out)   
!    Next, the variables used in debugging are declared. Including some variables used to identify the instant when the message is printed out
        EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
        INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout,                                    ! Output of erhandler (ANSYS own subroutine for displaying program errors) (inside)  
     &  ieqitr,                                    ! Current equilibrium iteration number (sub-in)                                    
     &  isubst,                                    ! Current substep number (sub-in)
     &  elId                                       ! Element number (sub-in)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!     Obtain from the old stiffness damage variable DkOld, the old energy damage variable DeOld
      ! If there was no damage in the previous substep (DkOld is 0), DeOld is also 0
      IF (DkOld < 1.0D-10) THEN
            DeOld = ZERO
      ! If there was complete degradation in the previous substep (DKold is 1), DeOld is also 1
      ELSE IF (DkOld > 0.9999999) THEN
            DeOld = ONE
      ! If there was damage in the previous substep, but not complete degradation (0 < DKold < 1), then calculate DeOld
      ELSE
            ! Firstly, obtain the one-dimensional opening displacement associated to the previous damage state. 
            ! Equation (14) of [2]
            deltaDOld = delta0Old*deltacOld / 
     & (deltacOld - DkOld*(deltacOld-delta0Old))
            ! Secondly, calculate the energy damage variable from the previous damage state. 
            ! Equation (20) of [2]
            DeOld = ONE - deltacOld*(1-DkOld)*KStiffness*deltaDOld
     & / (TWO*GcOld)
      END IF
      ! Generate an error message if the energy damage variable is not between 0 and 1
      IF (DeOld < ZERO .OR. DeOld > ONE) THEN
            errout = erhandler('USERELEM',5004,3,
     & 'ERROR in DeOld %F',(/DeOld/), ' ')
      END IF
!
!     Calculate the crack growth rate based on the energy release rate through the Paris' law
      ! Obtain the total energy release rate and mode mixity from the J-integral calculations
      G = GI + GII
      modemixity = GII/G
      ! Make the modemixity 0 in case it is negative. 
      ! This is done in case after observing sometimes a very small negative mode I contribution, which despite being close to 0 leads to errors because the mode mixity is raised to a fractional power for BKeta. Probably no longer necessary since now the J-integral contribution is made 0 if there is interpenetration
      IF (modemixity < ZERO) THEN
            modemixity = ZERO
      END IF
      ! Calculate the parameters of the Paris law for the current mode mixity. 
      ! Equations (44) from [2]      
      cParis = ((cParisII/(cParisI*cParism))**((modemixity)**2))* 
     & (cParism**(modemixity))* cParisI
      mParis = (modemixity)**2*(mParisII-mParisI-mParism) 
     & + modemixity*mParism + mParisI
      ! Obtain the global critical energy release rate through the BK criterion. 
      ! Equation (45) from [2]
      GcGlobal = GIc + (GIIc-GIc)*(modemixity)**BKeta
      ! Generate a warning message if the J-integral recorded is greater than the critical energy release rate
      IF (G > GcGlobal) THEN
            errout = erhandler('USERELEM',5005,2,
     & 'ERROR: G %F > GcGlobal %F
     & Substep %F, iteration %F and element %F.
     & Gs: %F, %F. Mode mixity: %F.
     & Properties: %F, %F, %F.',
     & (/G, GcGlobal,
     & isubst*ONE, ieqitr*ONE, elId*ONE,
     & GI, GII, modemixity, GIc, GIIc, BKeta/), ' ')
            ! This has been observed to occur in some specimens with unstable fatigue propagation. Full damage is then imposed, skipping the calculation of DkNew and moving directly to that of the dissipated work
            DkNew = ONE
            GO TO 100
      END IF
      ! Use the Paris' law to obtain the crack growth rate. 
      ! Equation (43) from [2], but without the threshold Gth and using G (which is Gmax since the maixmum load/displacement is being applied) instead of DeltaG
      dadN = cParis*(G/GcGlobal)**mParis
!
!     If there is damage. 
      IF ((delta-delta0) > 1.0D-10) THEN
!
!     Calculate the derivative of the energy damage variable with respect to the number of cycles (dDedNNew). 
      ! Equation (22) from [2], where the derivative with respect to the crack growth (dDeda) has been obtained beforehand through subroutine Rates
            dDedNNew = dDeda*dadN
! 
!     As explained in [2], the use of coarse meshes would lead to an unintended quasi-static damage. For this reason, a predictor-corrector approach is used
      ! In the predictor step, the unintended quasi-static damage is omitted and the energy damage variable is obtained by approximating the integral of the fatigue
      ! damage rate though the trapezoidal rule. Equation (30) in [2].
      ! First, obtain the derivative of the energy damage variable with respect to the number of cycles using the dDeda from the previous substep
            dDedNOld = dDedaOld*dadN
            ! Avoid undesired values of these derivatives
            IF (dDedNNew < ZERO) THEN
                  dDedNNew = ZERO
            END IF
            IF (dDedNOld < ZERO) THEN
                  dDedNOld = ZERO
            END IF
      ! Second, account for the corrected local time in the previous substep. 
            ! Equation (33) in [2]
            IF (Ncor > ZERO .AND. Ncor+Ns > Nloc) THEN
                  DeltaNcor = DeltaN * Nloc/(Ncor+Ns)
            ELSE
                  DeltaNcor = DeltaN
            END IF
      ! Finally, the energy damage variable is obtained. 
            ! Equation (30) substituted in (27)
            DeNew = DeOld + HALF*(dDedNOld + dDedNNew)*DeltaNcor   
            ! Avoid an energy damage variable greater than 1.0
            IF (DeNew > ONE) THEN
                  DeNew = ONE
            END IF
            ! Generate an error message if the new energy damage variable is lower than the previous one
            IF (DeNew < DeOld) THEN
                  errout = erhandler('USERELEM',5006,3,
     & 'ERROR: dDedNNew %F < dDedNOld %F
     & Substep %F, iteration %F and element %F.
     & Calculations: %F, %F, %F, %F, %F, %F.',
     & (/DeNew, DeOld,
     & isubst*ONE, ieqitr*ONE, elId*ONE,
     & dDedNOld, dDedNNew, DeltaNcor,
     & dDeda, dDedaOld, dadN/), ' ')
            END IF
!
!     In the corrector step, if unintended quasi-static damage occurs, such damage is calculated to later obtain the corrected local time
      ! First, the previously calculated energy damage due solely to fatigue is saved in a new variable
            DefatigueNew = DeNew
      ! Second, the threshold damage is calculated from the energy damage variable 
            ! Substituting equation (20) in (14) from [2]
            deltath = deltac + TWO*Gc*(DeNew-1)*(deltac-delta0)
     & /(deltac*KStiffness*delta0)
            ! This damage threshold is used to detected if non-intended quasi-static damage has occured
            ! If non-intended quasi-static damage has occured
            IF (delta > deltath) THEN
                  ! First, the stiffness damage variable is calculated as done under static loading
                  ! Equation (13) from [2]
                  DkNew = deltac*(delta - delta0) / 
     & (delta*(deltac - delta0))
                  ! Equation (12) from [2]
                  IF (DkNew > 0.9999999D0) THEN
                        DkNew = ONE
                  END IF
                  ! Second, the energy damage variable is obtained from the stiffness damage variable. 
                  ! Substituting equation (20) in (14) from [2]
                  DeNew = ONE - deltac*(ONE-DkNew)*KStiffness*
     & delta0*deltac/(TWO*Gc*(deltac-DkNew*(deltac-delta0)))
                  ! Avoid the energy damage variable being greater than 1.0
                  IF (DeNew > 0.9999999D0) THEN
                        DeNew = ONE
                  END IF
                  ! Third, the unintended evolution of quasi-static damage can then be calculated as the difference between both energy damage variables
                  DeltaDes = DeNew - DefatigueNew
                  ! Finally, the time variables can be updated
                  ! Equation (32) from [2]
                  IF ((dDedNOld+dDedNNew) == ZERO) THEN
                        Ns = Ns
                  ELSE
                        Ns = Ns + TWO*DeltaDes/
     & (dDedNOld+dDedNNew)
                  END IF
                  ! Equation (31) from [2]
                  Nloc = Nloc + DeltaN
                  ! Equation (33) from [2]
                  Ncor = Ncor + DeltaNcor
            ! If non-intended quasi-static damage has not occured
            ELSE
                  ! First, calculate the stiffness damage variable from the energy damage variable
                  ! Combination of equations (14) and (20) from [2]
                  DkNew = ((ONE-DeNew)*TWO*Gc*deltac 
     & - deltac**2*KStiffness*delta0) /
     & ((ONE-DeNew)*TWO*Gc*(deltac-delta0) 
     & - deltac**2*KStiffness*delta0)
                  ! Finally, the time variables are updated. Unlike if non-intended quasi-static damage occurs, Ns does not change
                  Ns = Ns 
                  Nloc = Nloc + DeltaN
                  Ncor = Ncor + DeltaNcor
            END IF 
!
!     It might be the case that DeltaN is very small for this specific crack, such that DeNew is very similar to DeOld
      ! This can cause that when calculating DkNew, the difference in Gc and deltac compared to the previous substep may lead to DkNew being slightly below DkOld, mainly due to changes in the mode mixity
      ! This is corrected next
            IF (DkNew < DkOld) THEN
                  ! If DeNew and DeOld are very similar, avoid healing by imposing DkNew to be the same as DkOld
                  IF ((DeNew-DeOld) < 1.0D-6) THEN
                        DkNew = DkOld
                  ! If DeNew is not so similar to DeOld, an error may exist in the model or the subroutine. But it can also be caused by highly unstable crack growth
                  ELSE
                        ! Generate an error message. Consider whether to leave this as a warning message and make DkNew = DkOld
                        errout = erhandler('USERELEM',5007,3,
     & 'ERROR: DkNew %F < DkOld %F
     & Substep %F, iteration %F and element %F
     & Inner variables %F, %F, %F, %F, %F, %F.
     & Number of cycles %F
     & De %F, %F, %F
     & Dk calculation %F, %F, %F, %F',
     & (/DkNew, DkOld,
     & isubst*ONE, ieqitr*ONE, elId*ONE,
     & delta, deltath, delta0, deltac, GI, (DkNew-DkOld)*1E6,
     & DeltaN,
     & DeOld, Denew, (DeNew-DeOld)*1E6,
     & Gc, deltac, GcOld, deltacOld/), ' ')
                  END IF 
            END IF
!
!     Detect full degradation and avoid the damage variables becoming greater than 1.0
            IF (delta > deltac .OR. DkNew > ONE) THEN
                  DkNew = ONE
                  DeNew = ONE
            END IF
!
!     If there is no degradation, the damage state does not change
      ! This ELSE comes from IF ((delta-delta0) > 1.0D-10) THEN further up
      ELSE
            DkNew = DkOld
            DeNew = DeOld
            DefatigueNew = DeNew
      END IF
!
!     Generate an error message if the damage variables are negative
      IF (DkNew < ZERO .OR. DeNew < ZERO) THEN
            errout = erhandler('USERELEM',5008,3,
     & 'ERROR: DkNew %F or DeNew %F
     & Substep %F, iteration %F and element %F.
     & Old: %F, %F
     & Is there damage? %F, %F',
     & (/DkNew, DeNew,
     & isubst*ONE, ieqitr*ONE, elId*ONE,
     & DkOld, DeOld,
     & delta, delta0/), ' ')
      END IF
!
!     Calculate the total dissipated specific work by solving equation (17) from [2] through the trapezoidal rule. It can also be seen as the area in figure 4 of [2]
      ! If there has been no damage
100   IF (DkNew < 0.1d-15) THEN
            ! In this case, the total dissipated specific work is a triangle of base delta and height KStiffness*delta
            wtot = HALF*KStiffness*delta**2
      ! If there has been damage
      ELSE
            ! This requires the calculation of the new and the old threshold damage variables and tractions using DeNew
            ! First, calculate the new values
            ! From a combination of equations (8) and (20) from [2]
            tauNew = TWO*(ONE-DeNew)*Gc/deltac
            ! Differentiate in case DeNew is still 0
            IF (DeNew > 1.0D-10) THEN
                  ! From a combination of equations (14) and (20) from [2]
                  deltathNew = TWO*Gc*(deltac*(DeNew-ONE) + 
     & delta0*(ONE-DeNew)) / (deltac*delta0*KStiffness) + deltac
            ELSE
                  deltathNew = delta0
            END IF
            ! Next, calculate the old values
            ! Differentiate between damage already having started or not
            IF (DkOld > 1.0D-10) THEN
                  ! From a combination of equations (14) and (20) from [2]
                  deltathOld = TWO*GcOld*(deltacOld*(DeOld-ONE) + 
     & delta0Old*(ONE-DeOld)) / (deltacOld*delta0Old*KStiffness) 
     & + deltacOld 
                  ! From a combination of equations (8) and (20) from [2]
                  tauOld = TWO*(ONE-DeOld)*GcOld/deltacOld
            ! If damage hasn't started in the previous step, then use the onset values as the previous values
            ELSE
                  wtotOld = HALF*tau0*delta0
                  tauOld = tau0
                  deltathOld = delta0
            END IF
            ! Finally, obtain the total dissipated specific energy by solving equation (17) from [2] via the trapezoidal rule
            wtot = wtotOld + HALF*(tauOld + tauNew) * 
     & MAX(ZERO, deltathNew - deltathOld)
      END IF
      ! Avoid it being negative
      IF (wtot < ZERO) THEN
            wtot = ZERO
      END IF
!
!     Use the calculated total dissipated specific work in the criterion for full damage
      ! Equation (35) from [2]
      ! Evaluate only if there is damage and full damage has not been reached already
      IF (Nloc > ZERO .AND. DkNew < ONE) THEN
            IF (wtot > G*(ONE + (Ncor+Ns-Nloc)/Nloc)) THEN
                  DkNew = ONE
                  DeNew = ONE
            END IF
      END IF
!
!     Debugging
!       errout = erhandler('USERELEM',1014,2,
!      & 'Fatigue-loading damage. Printing some values:
!      & Substep %F, iteration %F, element %F
!      & Deltas: %F, %F, %F, %F, %F, %F, %F, %F.
!      & Calculation De: %F, %F, %F, %F, %F, %F, %F, %F,
!      & %F, %F
!      & G: %F, %F, %F, %F, %F',
!      & (/ isubst*1.D0, ieqitr*ONE, elId*1.D0,
!      & deltac, delta0, delta, deltath, DkOld, DkNew,
!      & DeOld, DeNew,
!      & dDedNNew, dDedNOld, dDeda, dDedaOld, dadN,
!      & cParis*(G/GcGlobal)**mParis, G/GcGlobal, 
!      & DeltaN, cParis, mParis, 
!      & modemixity, GI, GII, G, GcGlobal/), 
!      & ' ')
!
      RETURN
      END SUBROUTINE DamageFatigue
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE InitiationVariables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Subroutine to calculate fatigue initiation variables
      ! These variables are: initiation slope, severity, initiation damage variable and cycle increment based on initiation
      SUBROUTINE InitiationVariables(sI, sII, mParisI, mParisII, deltaI,  
     & delta, deltaS, KStiffness, DkNew, sm, SeverityInit, tau0, 
     & DiNew, DiOld, DeltaN, DeltaNi,
     & ieqitr, isubst, elId)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      DOUBLE PRECISION
     &  mParisI,                                   ! Parameter m in Paris law for mode I (sub-in)
     &  mParisII,                                  ! Parameter m in Paris law for mode II (sub-in)
     &  sI,                                        ! Mode I initiation S-N curve slope (inside)
     &  sII,                                       ! Mode II initiation S-N curve slope (inside)
     &  cosI,                                      ! Mode I cosine director (inside)
     &  deltaI,                                    ! Mode I displacement jump (sub-in)
     &  delta,                                     ! Equivalent one-dimensional opening displacement (sub-in)
     &  cosII,                                     ! Mode II cosine director (inside)
     &  deltaS,                                    ! Shear (modes II and III) displacement jump (sub-in)
     &  sm,                                        ! Mode mixity S-N curve slope (sub-out)
     &  tau,                                       ! Interface traction in the current cycle increment (inside)
     &  KStiffness,                                ! Penalty stiffness. Initial stiffness of undamaged cohesive element (sub-in)
     &  DkNew,                                     ! Stiffness damage variable in the current substep (sub-in)
     &  SeverityInit,                              ! Severity: ratio of current traction over static strength (sub-out)
     &  Ni,                                        ! Number of cycles until initiation (inside)
     &  tau0,                                      ! Mixed mode equivalent interface strength (sub-in)
     &  DiNew,                                     ! Initiation damage variable in the current cycle increment (sub-out)
     &  DiOld,                                     ! Initiation damage variable in the previous cycle increment (sub-in)
     &  DeltaN,                                    ! Increment of cycles used in the calculations (inside)
     &  DeltaNi                                    ! Increment of cycles based on initiation for next substep (sub-out)
!    Next, the variables used in debugging are declared
        EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
        INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout,                                    ! Output of erhandler (ANSYS own subroutine for displaying program errors) (inside)   
     &  ieqitr,                                    ! Current equilibrium iteration number (sub-in)                                    
     &  isubst,                                    ! Current substep number (sub-in)
     &  elId                                       ! Element number (sub-in)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!     Approximate the initiation S-N curve slope from the Paris law slopes. This is based on [5]
      ! This is done if the slopes are 0 because the user in the ANSYS input file did not provide a value
      IF (sI == ZERO .AND. sII == ZERO) THEN
            sI = ONE/(TWO*mParisI)
            sII = ONE/(TWO*mParisII)
      END IF
!
      ! Calculate the slope of the S-N curve based on the displacement-based mode mixity, using the so-called cosine directors
      ! Equations (4) and (5) from [6]
      cosI = deltaI/delta
      cosII = deltaS/delta
      ! Equation (7) from [3] 
      sm = sI*cosI**2 + sII*cosII**2
!
!     Calculate the number of cycles until initiation for this stress state
      ! First, obtain the current traction
      ! Equation (8) from [2]
      tau = delta*KStiffness*(ONE-DkNew)
      ! Calculate the severity as the ratio between the current traction and the onset strength
      SeverityInit = tau/tau0
      ! Obtain the number of cycles to initiation
      ! Equation (8) from [3]
      Ni = 10.D0**((1-tau/tau0)/sm)
!
!     Calculate the new value of the initiation damage variable
      ! Discrete version of equation (11) from [3], combined as well with equations (9) and (10). Or equation (22) from [7]
      DiNew = DiOld + DeltaN/Ni
!
!     Calculate the cycles increment to initiation based on the current value of the initiation damage variable
      DeltaNi = (ONE-DiNew)*Ni
!
!     If already damaged, indicate an artificially high number of cycles to initiation to avoid it being considered
      IF (DkNew > ONE-1.0D-10) THEN
            DeltaNi = 1.0D10
      END IF
!
!       ! Debugging
!       errout = erhandler('USERELEM',1015,2,
!      & 'Inside Initiation subroutine
!      & Substep %F, iteration %F and element %F.
!      & Initiation variables: %F, %F, %F, %F.
!      & Calculation: %F, %F, %F, %F, %F, %F',
!      & (/isubst*ONE, ieqitr*ONE, elId*ONE,
!      & Ni, Ni/1e6, tau, tau0,
!      & cosI, cosII, sI, sII, sm,
!      & delta/), ' ')
!
      RETURN
      END SUBROUTINE InitiationVariables
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE TangentStiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Subroutine to calculate the element tangent stiffness matrix
      SUBROUTINE TangentStiffness(nDim, NIntP, IntP, nUsrDof, 
     & Loading, ConstStiff, DisplJump, KH, eStiff, Bglobal, dA)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      INTEGER 
     &  nDim,                                      ! Number of dimensions of the problem (sub-in)
     &  NIntP,                                     ! Number of integration points (sub-in)
     &  IntP,                                      ! Integration point number (sub-in)
     &  nUsrDof,                                   ! Number of DOFs of this element (matrix and load vector size) (sub-in)
     &  Loading(NIntP**2),                         ! Vector of flags indicating the loading or not of the cohesive element (sub-in)
     &  i,                                         ! Used in loops (inside)
     &  j,                                         ! Used in loops (inside) 
     &  k                                          ! Used in loops (inside)
        DOUBLE PRECISION
     &  TangentStiff(NIntP**2, nDim, nDim),        ! Material tangent constitutive matrix (inside)
     &  ConstStiff(NIntP**2, nDim, nDim),          ! Constitutive stiffness matrix (sub-in)
     &  DisplJump(NIntP**2, nDim),                 ! Displacement jump at integration point (sub-in)
     &  KH(NIntP**2),                              ! Parameter used in the calculation of the constitutive tangent stiffness matrix (sub-in)
     &  eStiff(nUsrDof,nUsrDof),                   ! Element tangent stiffness matrix (sub-inout)
     &  Bglobal(nDim, nUsrDof),                    ! B-matrix in global coordinates (sub-in)
     &  IntermTangStiff(nDim, nUsrDof),            ! Intermediate matrix in the calculation of the element tangent stiffness matrix (inside)
     &  dA                                         ! Area scale factor for current integration point value (sub-in)
!    Next, the variables used in debugging are declared
        EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
        INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout                                     ! Output of erhandler (ANSYS own subroutine for displaying program errors) (inside)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!     Initialise the variables
      DO i = 1, nDim
        DO j = 1, nDim
            TangentStiff(IntP, i, j) = ZERO
        END DO
        DO k = 1, nUsrDof
            IntermTangStiff(i, k) = ZERO
        END DO
      END DO
!
!     Calculate the material tangent stiffness matrix 
      ! Equation (51) from [1] 
      ! If the cohesive element is not being loaded, the tangent constitutive matrix coincides with the constitutive one
      IF (Loading(IntP) == 0) THEN
        DO i = 1, nDim
            DO j = 1, nDim
                TangentStiff(IntP, i, j) = ConstStiff(IntP, i, j)
            END DO
        END DO
      ! If the element is being loaded, the tangent constitutive matrix is calculated based on KH and the displacement jumps
      ELSEIF (Loading(IntP) == 1) THEN
        DO i = 1, nDim
            DO j = 1, nDim
                TangentStiff(IntP, i, j) = ConstStiff(IntP, i, j) -
     & KH(IntP) * DisplJump(IntP,j) * DisplJump(IntP,i)                   
            END DO
        END DO
        ! If penetration is observed, some components are not modified by damage
        IF (DisplJump(IntP, 3) < ZERO) THEN 
            DO i = 1, nDim
                TangentStiff(IntP, i, 3) = ConstStiff(IntP, i, 3)
                TangentStiff(IntP, 3, i) = ConstStiff(IntP, 3, i)
            END DO
        END IF
      END IF
!     
!     Calculate the element tangent stiffness matrix
      ! Discrete version of equation (68) from [1], and separating the products in two
      DO i = 1, nDim
        DO j = 1, nUsrDof
            DO k = 1, nDim
                IntermTangStiff(i, j) = IntermTangStiff(i, j) +
     & TangentStiff(IntP, i, k) * Bglobal(k, j)
            END DO
        END DO
      END DO
      DO i = 1, nUsrDof
        DO j = 1, nUsrDof
            DO k = 1, nDim
                eStiff(i, j) = eStiff(i, j) +
     & Bglobal(k, i) * IntermTangStiff(k, j) * dA
            END DO
        END DO
      END DO
!
!     Debugging
!       errout = erhandler('USERELEM',1016,2,
!      & 'TangentStiffness. Printing some values   
!      & %F, %F.',
!      & (/Loading(IntP)*ONE, eStiff(1, 1)/), ' ')
!     End of the subroutine
!
      RETURN
      END SUBROUTINE TangentStiffness
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE InternalForce %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Subroutine to calculate the internal force vector in global coordinates, as well as the tractions
      SUBROUTINE InternalForce(nDim, Traction, NIntP, IntP,  
     & nUsrDof, fInt, ConstStiff, DisplJump, Bglobal, dA)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      INTEGER 
     &  nDim,                                      ! Number of dimensions of the problem (sub-in)
     &  NIntP,                                     ! Number of integration points (sub-in)
     &  IntP,                                      ! Integration point number (sub-in)
     &  nUsrDof,                                   ! Number of DOFs of this element (matrix and load vector size) (sub-in)
     &  i,                                         ! Used in loops (inside)
     &  j                                          ! Used in loops (inside)    
        DOUBLE PRECISION
     &  fInt(nUsrDof),                             ! Internal force vector (sub-inout)
     &  Traction(NIntP**2, nDim),                  ! Traction matrix (sub-out)
     &  ConstStiff(NIntP**2, nDim, nDim),          ! Constitutive stiffness matrix (sub-in)
     &  DisplJump(NIntP**2, nDim),                 ! Displacement jump at integration point (sub-in)
     &  Bglobal(nDim, nUsrDof),                    ! B-matrix in global coordinates (sub-in)
     &  dA                                         ! Area scale factor for current integration point value (sub-in)
!    Next, the variables used in debugging are declared
        EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
        INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout                                     ! Output of erhandler (ANSYS own subroutine for displaying program errors) (inside)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!     Initialise the variables
      DO i = 1, nDim
        Traction(IntP, i) = ZERO
      END DO
!
!     Calculate the traction matrix in local coordinates
      ! Equation (18) from [1], but having incorporated already the stiffness damage variable in the constitutive stiffness matrix
      DO i = 1, nDim
        Traction(IntP, i) = ConstStiff(IntP, i, i)*
     & DisplJump(IntP, i)
      END DO
!
!     Calculate the material tangential stiffness matrix. 
      ! Discrete version of equation (67) from [1]
      DO i = 1, nUsrDof/2
        DO j = 1, nDim
            fInt(i) = fInt(i) + 
     & Bglobal(j, i) * Traction(IntP, j) * dA 
        END DO
        fInt(i + nUsrDof/2) = - fInt(i)
      END DO
!
! !     Debugging
!       errout = erhandler('USERELEM',1017,2,
!      & 'InternalForce. Printing some values   
!      & %F, %F, %F, %F.',
!      & (/Traction(IntP, 3), fInt(3),
!      & ConstStiff(IntP, 3, 3), DisplJump(IntP, 3)/)
!      & , ' ')
!
!     End of the subroutine
      RETURN
      END SUBROUTINE InternalForce
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE JIntegralCalc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Subroutine to calculate the J-integral element contributions
      SUBROUTINE JIntegralCalc(nDim, nNodes, xCur, NIntP,
     & elJIntegralI, elJIntegralII, elLength, ddelta3da, 
     & ddelta1da, Traction, IntPWeights, elId, isubst) 
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      INTEGER ::
     &  nDim,                                      ! Number of dimensions of the problem (in)
     &  nNodes,                                    ! Number of nodes of the element (in)
     &  NIntP,                                     ! Number of integration points (sub-in)
     &  i,                                         ! Used in loops
     &  j                                          ! Used in loops
      DOUBLE PRECISION ::
     &  elJIntegralI,                              ! Mode I J-integral contribution from this element (sub-out)
     &  elJIntegralII,                             ! Mode II J-integral contribution from this element (sub-out)
     &  xCur(nDim,nNodes),                         ! Nodal coordinates in current configuration (in)
     &  elLength,                                  ! Approximation of the element length (sub-out) 
     &  ddelta3da,                                 ! Derivative of the opening displacement with respect to the crack growth a (sub-in)
     &  ddelta1da,                                 ! Derivative of the mode II displacement with respect to the crack growth a (sub-in)
     &  Traction(NIntP**2, nDim),                  ! Traction matrix (sub-in)
     &  IntPWeights(NIntP)                         ! Weights for the integration points (sub-in)
!    Next, the variables used in debugging are declared. Including some variables used to identify the instant when the message is printed out
      EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
      INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout,                                    ! Output of erhandler (ANSYS own subroutine for displaying program errors)                                 
     &  isubst,                                    ! Current substep number (sub-in)
     &  elId                                       ! Element number (sub-in)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
! First, zero the element J-integral contributions in modes I and II
      elJIntegralI = ZERO
      elJIntegralII = ZERO
!
! Second, approximate the element length as the distance at the mid-surface between nodes 1 and 2
      elLength = SQRT( ((xCur(1, 2)+xCur(1, 6))/2 
     & - (xCur(1, 1)+xCur(1, 5))/2)**2 
     & +  ((xCur(2, 2)+xCur(2, 6))/2 
     & - (xCur(2, 1)+xCur(2, 5))/2)**2)
!
! Last, obtain the J-integral element contributions
      ! Equations (41) and (42) [corrected]
      ! The derivatives of the displacement jumps with respect to the crack growth have been previously calculated through subroutine Slopes. In this calculation, 
      ! the derivatives with respect to the crack growth coordinate xa are to be used, which can be approximated under self-similar growth to the derivatives with 
      ! respect to the crack growth a, but with opposite signs. Since equations (41) and (42) of [2] are incorrect, because a minus is missing, the calculation can 
      ! be made using the derivatives wrt a and not adding the minus sign
      ! Sum the terms from each integration point along one line in the crack growth direction xa
      j = 1
      DO i = 1, NIntP**2, NIntP 
            elJIntegralI = elJIntegralI + 
     & 0.5*elLength*ddelta3da*Traction(i,3)*IntPWeights(j)
            elJIntegralII = elJIntegralII + 
     & 0.5*elLength*ddelta1da*Traction(i,1)*IntPWeights(j)
            j = j + 1
      END DO
!
! ! Debugging
!       errout = erhandler('USERELEM',1018,2,
!      & 'Debugging J-integral calculation.
!      & Substep %F and element %F.
!      & Element length: %F   
!      & JI: %F, %F, %F, %F. 
!      & JII: %F, %F, %F, %F.',
!      & (/isubst*ONE, elId*ONE, 
!      & elLength,
!      & elJIntegralI, ddelta3da, 
!      & Traction(1,3), Traction(4,3), 
!      & elJIntegralII, ddelta1da, 
!      & Traction(1,1), Traction(4,1)/), ' ')
!
      RETURN
      END SUBROUTINE JIntegralCalc
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBROUTINE CrackTipIdentification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Subroutine to identify if the current element is a crack tip element and whether it already belongs to a crack or not
      SUBROUTINE CrackTipIdentification(NIntP, DkNews, CrackNumbers,
     & CrackTip, elId, cParisICrack, cParisI, cParisIICrack, 
     & cParisII, cParismCrack, cParism, mParisICrack, mParisI,
     & mParisIICrack, mParisII, mParismCrack, mParism, GIcCrack, 
     & GIc, GIIcCrack, GIIc, BKetaCrack, BKeta,
     & isubst, ieqitr)
      IMPLICIT NONE
!     Declare the variables used inside the subroutine
      INTEGER
     &  NIntP,                                     ! Number of integration points (sub-in)     
     &  CrackNumbers(1e5),                         ! Variable containing the crack number associated to each element (sub-in)
     &  CrackTip(10, 2),                           ! Variable containing the element situated in the crack tip (i,1) and the direction of propagation (i,2) for each crack i (sub-inout)
     &  elId,                                      ! Element number (sub-in)
     &  i,                                         ! Used in loops (inside)
     &  j,                                         ! Used in loops (inside)
     &  k                                          ! Used in loops (inside)
      DOUBLE PRECISION
     &  DkNews(NIntP**2),                          ! Variable containing the stiffness damage variable in the current substep for all integration points (sub-in)
     &  cParisICrack(10),                          ! Parameter c in Paris law for mode I for this crack. Defined in case there are multiple cracks growing in various interfaces with different interface properties (sub-out)
     &  cParisIICrack(10),                         ! Parameter c in Paris law for mode II for this crack (sub-out)
     &  cParismCrack(10),                          ! Parameter c in Paris law for fitting of mode interaction for this crack (sub-out)
     &  mParisICrack(10),                          ! Parameter m in Paris law for mode I for this crack (sub-out)
     &  mParisIICrack(10),                         ! Parameter m in Paris law for mode II for this crack (sub-out)
     &  mParismCrack(10),                          ! Parameter m in Paris law for fitting of mode interaction for this crack (sub-out)
     &  GIcCrack(10),                              ! Mode I critical energy release rate / fracture toughness for this crack (sub-out)
     &  GIIcCrack(10),                             ! Mode II critical energy release rate / fracture toughness for this crack (sub-out)
     &  BKetaCrack(10),                            ! Benzeggagh–Kenane criterion coefficient for this crack (sub-out)
     &  cParisI,                                   ! Parameter c in Paris law for mode I (sub-in)
     &  cParisII,                                  ! Parameter c in Paris law for mode II (sub-in)
     &  cParism,                                   ! Parameter c in Paris law for fitting of mode interaction (sub-in)
     &  mParisI,                                   ! Parameter m in Paris law for mode I (sub-in)
     &  mParisII,                                  ! Parameter m in Paris law for mode II (sub-in)
     &  mParism,                                   ! Parameter m in Paris law for fitting of mode interaction (sub-in)
     &  GIc,                                       ! Mode I critical energy release rate / fracture toughness (sub-in)
     &  GIIc,                                      ! Mode II critical energy release rate / fracture toughness (sub-in)
     &  BKeta                                      ! Benzeggagh–Kenane criterion coefficient (sub-in)
!    Next, the variables used in debugging are declared. Including some variables used to identify the instant when the message is printed out
      EXTERNAL ::
     &  erhandler                                  ! ANSYS own subroutine for displaying program errors
      INTEGER ::
     &  erhandler,                                 ! ANSYS own subroutine for displaying program errors
     &  errout,                                    ! Output of erhandler (ANSYS own subroutine for displaying program errors)                                 
     &  isubst,                                    ! Current substep number (sub-in)
     &  ieqitr                                     ! Current equilibrium iteration number (in)
!     Declare some parameters
      DOUBLE PRECISION ::
     & ZERO, ONE, TWO, FOUR, HALF
      PARAMETER (ZERO = 0.d0,
     & ONE = 1.d0, 
     & TWO = 2.d0,
     & FOUR = 4.d0,
     & HALF = 0.5d0)
!
!     Loop through the integration points
      LoopCrackTip: DO i = 1, NIntP**2
            DO j = 1, NIntP**2
                  ! The element is considered to be at the crack tip if one integration point is damaged, but not fully, and another one is undamaged
                  ! If one integration point is damaged, but not fully
                  IF (DkNews(i) < ONE .AND. DkNews(i) > 1.0D-6) THEN
                        ! If another integration point is undamaged
                        IF (DkNews(j) < 1.0D-6) THEN
                              ! If the element does belong to a crack
                              IF (CrackNumbers(elId) .NE. 1) THEN
                                    ! Save the element number for this crack in the variable CrackTip 
                                    CrackTip(CrackNumbers(elId),1) = 
     & elId
                                    k = CrackNumbers(elId) 
                              ! If the element does not belong to a crack 
                              ELSE
                                    ! Loop over the possible crack numbers
                                    LoopCrackNum: DO k = 2, 
     & SIZE(CrackTip, DIM = 1)
                                          ! Find the crack number which has not been used yet
                                          IF (CrackTip(k,1) == 0) THEN
                                                ! Save the element number for unused crack number k
                                                CrackTip(k,1) = elId
                                                EXIT LoopCrackNum
                                          END IF
                                    END DO LoopCrackNum
                              END IF
                              !
                              ! Calculate the propagation direction based on in which side the damage variable is greater: 0 (lower coordinate), 1 (higher)
                              IF (DkNews(1) > DkNews(NIntP**2)) THEN
                                    CrackTip(k,2) = 1
                              ELSE
                                    CrackTip(k,2) = 0
                              END IF
                              !
                              ! After saving the crack number and the propagation direction, assign the interface properties of this crack based on those of its crack tip element
                              cParisICrack(k) = cParisI
                              cParisIICrack(k) = cParisII
                              cParismCrack(k) = cParism
                              mParisICrack(k) = mParisI
                              mParisIICrack(k) = mParisII
                              mParismCrack(k) = mParism
                              GIcCrack(k) = GIc
                              GIIcCrack(k) = GIIc
                              BKetaCrack(k) = BKeta
                              !
!                               ! Debugging
!                               errout = erhandler('USERELEM',1019,2,
!      & 'Debugging crack tip identification.
!      & Substep %F and element %F.
!      & NumCrack: %F, GIIc: %F.',
!      & (/isubst*ONE, elId*ONE, 
!      & k*ONE, GIIcCrack(k)/), ' ')
                              !
                              ! Exit the loop
                              EXIT LoopCrackTip
                        END IF
                  END IF
            END DO
      END DO LoopCrackTip
!
!
      RETURN
      END SUBROUTINE CrackTipIdentification
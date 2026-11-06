BIMOND homogeneous evolution
============================

The experimental ``camb.bimond`` module evolves the homogeneous equations of a
specified bimetric modified Newtonian dynamics (BIMOND) action. The ordinary
sector contains pressureless baryons and radiation. The twin matter sector is
empty. Both metrics remain dynamical. The module provides a background initial
value problem, not BIMOND CMB spectra or density and velocity transfer functions.
The standard ``CAMBparams`` and Einstein--Boltzmann equations remain separate.

Action and units
----------------

We select the action of `Milgrom, equation 15
<https://arxiv.org/html/0912.0790>`_ with :math:`\alpha=-1`,
:math:`\beta=1`, :math:`f(\kappa)=1`, and
:math:`Z=-g^{\mu\nu}\Upsilon_{\mu\nu}/(2a_0^2)`.
The Ricci convention in equation 8 of the source reverses the usual modern GR
sign. Removing the Einstein--Hilbert boundary terms gives a negative ordinary
scale-factor kinetic term and a positive twin scale-factor kinetic term.

The flat metrics and reduced fluid action are

.. math::

   ds_g^2=-N^2dt^2+a^2d\mathbf{x}^2,\qquad
   ds_{\hat g}^2=-X^2dt^2+b^2d\mathbf{x}^2,

   8\pi G L=-\frac{3a\dot a^2}{N}
      +\frac{3b\dot b^2}{X}
      +a_0^2\sqrt{NX}(ab)^{3/2}M(Z)
      -8\pi G N(\rho_{b0}+\rho_{r0}/a).

The dimensional expression uses :math:`c=1`. The implementation uses an explicit
reference rate :math:`H_{\rm ref}`, dimensionless time
:math:`\tau=H_{\rm ref}t`, and parameters

.. math::

   \bar a_0=\frac{a_{0,\rm physical}}{cH_{\rm ref}},\qquad
   B=\frac{8\pi G\rho_{b0}}{H_{\rm ref}^2},\qquad
   R_r=\frac{8\pi G\rho_{r0}}{H_{\rm ref}^2}.

The densities are mass densities, including radiation energy density divided by
:math:`c^2`. The reference rate is a choice of units. The reference rate need not
equal a fitted present-day Hubble rate. Baryon and radiation densities evolve as
:math:`\rho_b=\rho_{b0}a^{-3}` and :math:`\rho_r=\rho_{r0}a^{-4}`.

Free function
-------------

The positive-:math:`Z` branch uses the simple QUMOND interpolation function

.. math::

   \nu(y)=\frac12+\sqrt{\frac14+\frac1y}.

We specify the negative branch by the signed continuation

.. math::

   M(Z)=M_0+\operatorname{sgn}(Z)\int_0^{|Z|}
         [\nu(\sqrt x)-1]dx.

The continuation is an additional research assumption. Galactic QUMOND does not
select the negative branch or :math:`M_0`. The constant :math:`M_0` contributes
to both metric equations through the mixed volume element and cannot be
identified directly with an ordinary cosmological constant.

For :math:`Z\ne0`,

.. math::

   M'=\frac{\sqrt{1+4/\sqrt{|Z|}}-1}{2},\qquad
   M''=-\frac{\operatorname{sgn}(Z)}
      {2|Z|^{3/2}\sqrt{1+4/\sqrt{|Z|}}}.

The implementation evaluates the primitive analytically. A convergent small-
:math:`\sqrt{|Z|}` series avoids subtractive cancellation. The numerical domain
excludes :math:`|Z|\le Z_{\min}`. The exclusion terminates evolution near the
nonanalytic surface and never smooths or replaces the selected free function.

Equations and constraint
------------------------

Define :math:`n=\ln N`, :math:`q=(u,v,w)=(\ln a,\ln b,\ln(X/N))`, and
:math:`(U,V,W)=dq/d\tau`. The common lapse derivative cancels exactly. With

.. math::

   s=e^{2(v-u-w)},\quad k=\frac{3}{2\bar a_0^2},\quad
   Q=2U^2-(3+s)UV+(1+s)V^2+(s-1)VW,

   T=-3e^{3u}U^2+3e^{3v-w}V^2,\quad
   C=\bar a_0^2 e^{3(u+v)/2+w/2},\quad R=B+R_r e^{-u},

the dimensionless reduced Lagrangian is

.. math::

   \mathcal L=e^{-n}T+e^n C M(ke^{-2n}Q)-e^nR.

After varying the common lapse, set :math:`N=1`. The three remaining
Euler--Lagrange equations and the common-lapse constraint become

.. math::

   H_{ij}\ddot q_j=\mathcal L_{,q_i}
         -\mathcal L_{,\dot q_i q_j}\dot q_j,
   \qquad H_{ij}=\mathcal L_{,\dot q_i\dot q_j},

   E=T+C[2ZM'-M]+R=0,\qquad Z=kQ.

Conservation of an arbitrary nonzero :math:`E` is insufficient. Initial states
must satisfy :math:`E=0`. The numerical integrator evolves all three equations
without projection onto the constraint. The reported constraint residual
measures the error in the omitted lapse equation.

Initial conditions may specify :math:`q,U,V` and solve
:math:`2ZM'-M=-(T+R)/C` for a selected signed branch. The resulting
:math:`Z` fixes :math:`W` when :math:`V(s-1)\ne0`. The initial constraint
solver does not adjust the model densities or :math:`M_0`.

Numerical and scientific domain
-------------------------------

The current ODE system requires a nonsingular velocity Hessian. In particular,
:math:`s=1` removes the :math:`w` row and column even for nonzero :math:`Z`.
The interface rejects singular or ill-conditioned Hessians and states near
:math:`Z=0`. A singular branch requires a separate constrained analysis.
The signs of the homogeneous Hessian do not establish the presence or absence
of ghosts in spatially varying perturbations.

The domain checks also apply to internal Runge--Kutta trial stages. A failed
trial stage can precede an accepted trajectory endpoint. Repeating the same
interval with a smaller ``max_step`` distinguishes overshoot from a persistent
approach to a singular surface. A failed trial alone does not locate a physical
singularity.

An unequal-metric background also requires a derivation of the local weak-field
limit before coupling to the existing lagRamses QUMOND operator. A regular
background integration alone does not validate the nonrelativistic matching.
Linear scalar perturbations, regular initial modes, relativistic stability,
recombination coupled to the selected expansion, and observable transfer
functions remain implementation requirements.

The source-level checks compare the Fortran action derivatives with independent
symbolic differentiation, verify the signed primitive by quadrature, and evolve
constraint-satisfying initial states over the same interval at increasing
accuracy. Such checks establish the numerical implementation within the stated
domain, not observational viability.

Tensor perturbations of the selected action
-------------------------------------------

The principal tensor action supplies a separate physical stability test.
For one transverse traceless polarization with
:math:`e_{ij}e_{ij}=2`, define :math:`\mathbf H=(h,\hat h)` and

.. math::

   r=b/a,\qquad s=r^2/X^2,\qquad t=r^3/X,\qquad
   m=\sqrt X\,r^{3/2}M'(Z)>0.

Direct contraction of the connection differences gives

.. math::

   S^{(2)}_{\rm principal}=\frac1{32\pi G}\int dt\,d^3x\,a^3
      \left[\dot{\mathbf H}^{\mathsf T}K\dot{\mathbf H}
      -a^{-2}\partial_i\mathbf H^{\mathsf T}G_T\partial_i\mathbf H\right],

   K=\begin{pmatrix}1-m&ms\\ms&-t-m(2s-1)\end{pmatrix},\qquad
   G_T=\begin{pmatrix}1-m&m\\m&-t/s-m\end{pmatrix}.

The first-order tensor variation of :math:`Z` vanishes by isotropy and
tracelessness. The coefficient :math:`M''` therefore contributes no principal
quadratic tensor term. The independent verification script constructs the
connections from metric derivatives and compares the contracted action with
both matrices.

Positive tensor kinetic energy requires :math:`K>0`, hence :math:`m<1`.
However,

.. math::

   \det G_T=-\frac ts(1-m)-m<0

for positive :math:`s,t,m` and :math:`m<1`. Such a background has a tensor
gradient instability. A nonsingular kinetic matrix outside the positive-definite
region instead has a tensor mode with negative kinetic energy relative to
ordinary matter. A singular kinetic matrix requires a separate constrained
analysis. No regular region of the selected action and continuation supplies
two propagating tensor modes free of both instabilities.

The obstruction applies to both signs of :math:`Z`, because the specified
continuation has :math:`M'>0` on both branches. The obstruction does not exclude
other BIMOND actions. The implemented homogeneous module remains a research
realization of the declared action and cannot supply a demonstrated viable
BIMOND replacement for the standard cosmology. A different action requires a
new derivation of the background, perturbations, and local force limit.

API
---

An arbitrary dimensionless example constructs a regular local initial state
without changing the supplied model densities.

.. code-block:: python

   import numpy as np
   from camb.bimond import (
       BIMONDParameters, constraint_initial_velocity, integrate_background,
       tensor_principal_matrices,
   )

   parameters = BIMONDParameters(a0bar=1.0, baryons=0.3, radiation=0.02)
   q = np.array([0.0, 0.2, -0.1])
   velocity = constraint_initial_velocity(
       q, 1.0, 0.65, parameters, branch=1, z_bracket=(1e-8, 1e6)
   )
   solution = integrate_background(
       parameters, np.concatenate((q, velocity)), (0.0, 0.01),
       sample_times=np.linspace(0.0, 0.01, 21), rtol=1e-11, atol=1e-13,
       max_step=0.01,
   )
   tensors = tensor_principal_matrices(q, solution.z[0], parameters)

The example parameters are numerical fixtures, not a cosmological fit.

Validation on 7 November 2026
----------------------------

Grammar Slurm job ``1274024`` built the complete lagCAMB library with gfortran
13.2.0. Ten BIMOND tests, nine independent symbolic tensor identities, and all
18 default CAMB regression tests passed. The source archive starts from lagCAMB
``d65514d5b99390d5f96cdb8546e9a45e7b422ae9`` and forutils
``1e4a56713c1b42c3d007fa175b8a32ccf7803cb3``. The archive SHA256 is
``8a1b5249e2abaa4647c6b688acd3ac973bbb61d656d783263bd9a2147cd344a4``.

The maximum scaled difference between the Fortran derivatives and independent
symbolic differentiation was :math:`1.89\times10^{-15}`. Independent quadrature
of the primitive agreed to :math:`6.67\times10^{-16}` relative error over
:math:`10^{-10}\le|Z|\le10^{10}`. On the fixed interval
:math:`0\le\tau\le0.01`, the tight integrations had maximum normalized
constraint errors of :math:`2.18\times10^{-16}` for the regular positive branch
and :math:`2.85\times10^{-14}` for the negative branch. Decreasing tolerances
reduced the negative-branch errors. The positive-branch runs used identical
steps and reached roundoff precision, so the positive fixture does not measure
a convergence order.

A separate positive-branch initial state with :math:`V=0.4` reached the
homogeneous Hessian exclusion near :math:`\tau=0.00320474627`. Grammar job
``1274015`` reproduced the approach using maximum steps of :math:`10^{-3}` and
:math:`10^{-4}`. The failed trajectory remains an explicit rejection test.
The successful positive fixture uses :math:`V=0.65` over the original interval.

Both successful numerical fixtures have a negative tensor kinetic eigenvalue.
The source-independent connection contraction confirms the broader tensor
obstruction for the specified action. Numerical integration success does not
remove the physical instability.

Compact reports reside in ``fortran/tests/bimond_validation_2026-11-07/``.

.. automodule:: camb.bimond
   :members:

Scalar and vector variational equations
--------------------------------------

The separate ``camb.bimond_symbolic`` module generates the quadratic action for
general Einstein--Hilbert coefficients :math:`\alpha,\beta` and :math:`f=1`.
The generator expands both inverse metrics, both connections, the ordinary
contraction of :math:`\Upsilon_{\mu\nu}`, and the three volume factors through
second order. The resulting action retains background terms and terms with one
derivative as well as the principal terms. No gauge or perturbation constraint
is imposed during generation.

The scalar coordinates are
:math:`(\phi,S,\zeta,\ell,\hat\phi,T,\hat\zeta,\hat\ell)`.
The ordinary lapse is :math:`e^\phi`, the twin lapse is :math:`Xe^{\hat\phi}`,
and the ordinary spatial metric is
:math:`a^2\operatorname{diag}(e^{2\zeta},e^{2\zeta},e^{2\ell})`.
The twin spatial metric uses :math:`b,\hat\zeta,\hat\ell`.
The contravariant longitudinal shifts are :math:`S,T`. Both shifts and both
lapse perturbations remain in the generated action.

The vector coordinates are :math:`(P,S,Q,T)` for a wave along the third spatial
axis and one transverse component along the first axis. The off-diagonal
spatial metric perturbations are :math:`a^2P` and :math:`b^2Q`, with
:math:`P=\partial_z F_x` and :math:`Q=\partial_z\hat F_x`. The contravariant
transverse shifts are :math:`S,T`. Quadratic ADM shift terms in the temporal
metric components are retained.

Grammar job ``1274098`` verified 19 identities for the gravitational generator.
The checks recover the general-lapse homogeneous action, the common lapse
constraint direction, both scalar principal gauge directions, and the
independently derived vector action. The full vector action has no shift
velocity Hessian and no antisymmetric shift symplectic pair. The scalar output
contains eight-field principal matrices and the complete unreduced quadratic
action. Eigenvalues of the unreduced scalar Hessian do not constitute a scalar
stability test.

Positive Einstein--Hilbert candidate
-----------------------------------

A subsequent research candidate uses :math:`\alpha=\beta=1`, conventional
positive matter actions, and

.. math::

   M'(Z)=\frac2{4+\sqrt{|Z|}},\qquad
   M(Z)=M_0+\operatorname{sgn}(Z)
       \left[4\sqrt{|Z|}-16\ln\left(1+\frac{\sqrt{|Z|}}4\right)\right].

The interaction scalar still contracts with the ordinary inverse metric alone.
Equal Einstein--Hilbert coefficients therefore do not make the action
exchange-symmetric. The fixed numerical ``camb.bimond`` background model
retains the earlier :math:`\alpha=-1,\beta=1` action. The symbolic generator
keeps the two candidates distinct.

For the positive-coefficient candidate, let
:math:`u=\dot F-S`, :math:`v=\dot{\hat F}-T`, and
:math:`d=\dot F-\dot{\hat F}`. The principal vector action for one transverse
component is

.. math::

   L_{V,\mathrm{principal}}=\frac{a^3k^2}{32\pi G}
       [(1-m)u^2+2muv+(t-m)v^2-2m(1-s)dv].

Eliminating the principal shift constraints gives

.. math::

   L_{V,\mathrm{reduced}}=\frac{a^3k^2}{32\pi G}C_Vd^2,\qquad
   C_V=-\frac{m^2(1-s)^2(1-m)}{t-m(1+t)}.

The near-equal metric-coefficient region has :math:`t-m(1+t)>0` for the
specified function and finite nonzero :math:`Z`. A generic displacement with
:math:`s\ne1` therefore has a negative relative-vector kinetic coefficient.
The region

.. math::

   0<s<1,\qquad 0<m<1,\qquad ms<t(1-m)<m

can pass the principal tensor and vector kinetic conditions. The interaction
and background equations determine :math:`m,s,t` jointly. The inequalities do
not establish a cosmological solution. The vector has zero principal spatial
speed, so lower-order terms remain necessary for finite-frequency evolution.
The surfaces :math:`s=1` and :math:`t-m(1+t)=0` require a separate constrained
analysis.

The optional scalar matter completion supplies irrotational dust and a
radiation scalar,

.. math::

   L_d=-\frac12\sqrt{-g}\,\lambda
       (g^{\mu\nu}\partial_\mu T_d\partial_\nu T_d+1),\qquad
   L_r=\sqrt{-g}\,\kappa Y^2,\quad
   Y=-\frac12g^{\mu\nu}\partial_\mu\chi\partial_\nu\chi.

Positive dust density and :math:`\kappa>0,Y>0` provide healthy material
kinetics for a constrained scalar calculation. The radiation scalar has
:math:`p_r=\rho_r/3` and sound speed squared :math:`1/3`. The matter completion
serves a kinetic and constraint calculation. Photon and neutrino hierarchies,
anisotropic stress, baryon thermodynamics, and recombination remain necessary
for physical transfer functions.

Constrained scalar result at two local backgrounds
-------------------------------------------------

Grammar Slurm job ``1274137`` repeated all 19 gravitational identities after
adding the matter action and evaluated two local backgrounds at 80-digit
precision. The backgrounds use :math:`M_0=0`, :math:`\bar a_0=1`,
:math:`a=1`, :math:`b=1/\sqrt2`, :math:`X=\sqrt2`, :math:`U=1`, and
:math:`V=1/2`. The selected values
:math:`Z=\pm[4(\sqrt2-1)]^2` give :math:`s=t=m=1/4`, inside the principal
tensor/vector-admissible region. Radiation density is :math:`0.02`, with
:math:`\dot{\bar\chi}=1` and :math:`\kappa=2/75` in the declared units.

The homogeneous constraint fixes positive dust densities of approximately
:math:`2.56018` and :math:`3.77482` for positive and negative :math:`Z`.
The homogeneous equations fix both scale-factor accelerations and the relative
lapse acceleration. The full matter-completed linear action vanishes modulo
boundaries on both local solutions. The maximum normalized tadpole residuals
are :math:`4.01\times10^{-80}` and :math:`3.61\times10^{-79}`.

The dust clock fixes the common time coordinate, and :math:`\ell=\zeta` fixes
the common scalar spatial coordinate. Variation of the dust multiplier imposes
:math:`\phi=0`. The ordinary lapse equation reconstructs the dust density
perturbation. Both shifts pass the auxiliary-variable checks. Integration by
parts retains background time derivatives before eliminating the shifts at
finite wavenumber. The remaining variables are
:math:`(\zeta,\hat\phi,\hat\zeta,\hat\ell,\delta\chi)`.

At dimensionless wavenumbers :math:`k=10,100,1000`, both backgrounds have a
nonsingular five-variable kinetic matrix with two negative eigenvalues.
At :math:`k=100`, the negative eigenvalues are approximately
:math:`(-64936.24,-0.111006)` for positive :math:`Z` and
:math:`(-27576.50,-0.708438)` for negative :math:`Z`.
The eigenvalue magnitudes depend on coordinate normalization and are not
instability rates. The number of negative physical kinetic directions is
invariant under nonsingular changes of field variables.

The radiation kinetic block is exactly :math:`K_{\chi\chi}=0.04` with zero
kinetic cross terms. The saved negative-direction certificates have zero
radiation component. The two negative directions therefore occupy the metric
and dust-clock sector. The nonsingular reduced kinetic matrix leaves no
additional local primary constraint capable of removing the negative
directions. Both tested backgrounds contain scalar ghosts. The result excludes
the two local backgrounds, not the whole positive-coefficient action or other
BIMOND theories. No characteristic growth rate or transfer function has been
computed.

The complete record is
``fortran/tests/bimond_validation_2026-11-07/scalar_constraints.json``.
The generator and standalone audit can be run on an allocated compute node
with

.. code-block:: bash

   python fortran/tests/verify_bimond_scalar_vector.py camb/bimond_symbolic.py NEW_RESULT_DIRECTORY
   python fortran/tests/audit_bimond_scalar_constraints.py --output NEW_SCALAR_REPORT.json

Both commands preserve existing result files. The audit rejects singular
inverses and reports unresolved null sectors without a stability verdict.

.. automodule:: camb.bimond_symbolic
   :members: QuadraticAction, quadratic_action, principal_part

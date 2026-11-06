! Homogeneous alpha=-1, beta=1 BIMOND, f(kappa)=1, no twin matter.
! See docs/source/bimond.rst for the action, units, and unresolved perturbations.
module BimondBackground
    use iso_c_binding, only: c_double, c_int
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: bimond_function, bimond_background

contains

    subroutine bimond_function(z, m0, zmin, values, status) bind(C, name="camb_bimond_function")
        real(c_double), intent(in) :: z, m0, zmin
        real(c_double), intent(out) :: values(3)
        integer(c_int), intent(out) :: status
        real(c_double) :: x, y, t, r, primitive, coefficient, series
        integer :: n

        values = 0
        status = 1
        if (.not.all(ieee_is_finite([z, m0, zmin]))) return
        if (zmin <= 0) return
        status = 2
        if (abs(z) <= zmin) return
        x = abs(z)
        y = sqrt(x)
        t = sqrt(y)
        r = sqrt(1 + 4/y)
        ! Integrate M' from zero. The series avoids cancellation at small |Z|.
        if (y < 0.1_c_double) then
            coefficient = 1
            series = 4._c_double/3
            do n = 1, 12
                coefficient = coefficient*(1.5_c_double-n)*y/(4*n)
                series = series + 4*coefficient/(3+2*n)
            end do
            primitive = t*y*series - x/2
        else
            primitive = y*(r + 2/(r+1)) - 4*asinh(t/2)
        end if
        values(1) = m0 + sign(primitive, z)
        values(2) = 2/(y*(r+1))
        values(3) = -sign(1._c_double, z)/(2*r*y*x)
        status = 1
        if (.not.all(ieee_is_finite(values))) return
        status = 0
    end subroutine bimond_function

    subroutine bimond_background(q, velocity, parameters, zmin, scalars, gradient, mixed, hessian, status) &
            bind(C, name="camb_bimond_background")
        ! q=(ln a, ln b, ln(X/N)); velocity=dq/d(H_ref t), evaluated at N=1.
        ! parameters=(a0/(c H_ref), 8piG rho_b0/H_ref^2, 8piG rho_r0/H_ref^2, M0).
        ! scalars=(L, E, constraint normalization, Z).
        ! mixed(i,j)=d^2 L/dvelocity_i dq_j; hessian(i,j)=d^2 L/dvelocity_i dvelocity_j.
        real(c_double), intent(in) :: q(3), velocity(3), parameters(4), zmin
        real(c_double), intent(out) :: scalars(4), gradient(3), mixed(3,3), hessian(3,3)
        integer(c_int), intent(out) :: status
        real(c_double) :: a3, b3x, s, k, amplitude, radiation, density, kinetic, z
        real(c_double) :: u, v, w, functions(3), direction(3), cq(3), tq(3), zq(3), zv(3)
        real(c_double) :: zvv(3,3), zvq(3,3), tvq(3,3), tvv(3,3), ds(3), dzvds(3)
        integer :: i, j

        scalars = 0
        gradient = 0
        mixed = 0
        hessian = 0
        status = 1
        if (.not.all(ieee_is_finite(q))) return
        if (.not.all(ieee_is_finite(velocity))) return
        if (.not.all(ieee_is_finite(parameters))) return
        if (parameters(1) <= 0 .or. any(parameters(2:3) < 0)) return
        ! Avoid exponential overflow before evaluating the action.
        if (maxval(abs(q)) > 100._c_double .or. maxval(abs(velocity)) > 1.e50_c_double) return
        a3 = exp(3*q(1))
        b3x = exp(3*q(2)-q(3))
        s = exp(2*(q(2)-q(1)-q(3)))
        k = 1.5_c_double/parameters(1)**2
        amplitude = parameters(1)**2*exp(1.5_c_double*(q(1)+q(2))+q(3)/2)
        radiation = parameters(3)*exp(-q(1))
        density = parameters(2) + radiation
        u = velocity(1)
        v = velocity(2)
        w = velocity(3)
        kinetic = -3*a3*u*u + 3*b3x*v*v
        z = k*(2*u*u-(3+s)*u*v+(1+s)*v*v+(s-1)*v*w)
        call bimond_function(z, parameters(4), zmin, functions, status)
        if (status /= 0) return
        scalars(1) = kinetic + amplitude*functions(1) - density
        scalars(2) = kinetic + amplitude*(2*z*functions(2)-functions(1)) + density
        scalars(3) = 3*a3*u*u + 3*b3x*v*v + &
            amplitude*(abs(2*z*functions(2))+abs(functions(1))) + density
        scalars(4) = z

        direction = [-1._c_double, 1._c_double, -1._c_double]
        ds = 2*s*direction
        cq = amplitude*[1.5_c_double, 1.5_c_double, 0.5_c_double]
        tq = [-9*a3*u*u, 9*b3x*v*v, -3*b3x*v*v]
        zq = k*v*(v-u+w)*ds
        zv = k*[4*u-(3+s)*v, -(3+s)*u+2*(1+s)*v+(s-1)*w, (s-1)*v]
        dzvds = k*[-v, -u+2*v+w, v]
        gradient = tq + cq*functions(1) + amplitude*functions(2)*zq
        gradient(1) = gradient(1) + radiation
        tvq = 0
        tvq(1,1) = -18*a3*u
        tvq(2,2) = 18*b3x*v
        tvq(2,3) = -6*b3x*v
        tvv = 0
        tvv(1,1) = -6*a3
        tvv(2,2) = 6*b3x
        zvv = 0
        zvv(1,1) = 4*k
        zvv(1,2) = -(3+s)*k
        zvv(2,1) = zvv(1,2)
        zvv(2,2) = 2*(1+s)*k
        zvv(2,3) = (s-1)*k
        zvv(3,2) = zvv(2,3)
        do j = 1, 3
            do i = 1, 3
                zvq(i,j) = dzvds(i)*ds(j)
                mixed(i,j) = tvq(i,j) + zv(i)*cq(j)*functions(2) + &
                    amplitude*(functions(3)*zv(i)*zq(j)+functions(2)*zvq(i,j))
                hessian(i,j) = tvv(i,j) + &
                    amplitude*(functions(3)*zv(i)*zv(j)+functions(2)*zvv(i,j))
            end do
        end do
        status = 1
        if (.not.all(ieee_is_finite(scalars))) return
        if (.not.all(ieee_is_finite(gradient))) return
        if (.not.all(ieee_is_finite(mixed))) return
        if (.not.all(ieee_is_finite(hessian))) return
        status = 0
    end subroutine bimond_background
end module BimondBackground

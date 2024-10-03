module emission
	implicit none
	integer, parameter :: dp = kind(8)
	real(dp), parameter :: c_cs = 1.199985_dp ! shottky const
    real(dp), parameter :: c_qe  = 1.6d-19
    real(dp), parameter :: c_me  = 9.11d-31
    real(dp), parameter :: c_hp  = 1.054d-34
    real(dp), parameter :: c_zs = 1.6118d14
    real(dp), parameter :: c_kb = 1.0_dp/11600.0_dp

contains
	function j_c(ts, f, ef, w)
        real(dp) :: ts, f, ef, w
        real(dp) :: j_c

        j_c = j_rds(ts, f, ef, w)+j_fn(ts, f, ef, w) 

    end function j_c

    function j_rds(ts, f, ef, w)
        real(dp) :: ts, f, ef, w
        real(dp) :: j_rds
		
		j_rds = 1.2d6*ts**2 * exp( -(w - c_cs*sqrt(f)) / (c_kb*ts) )

    end function j_rds

    function j_fn(ts, f, ef, w)
        real(dp) :: ts, f, ef, w
        real(dp) :: j_fn
        
        real(dp) :: a, b, y, v, t
        
        a = 1.54d12
        b = 6.83_dp
        y = c_cs*sqrt(f)/w
        v = 1.0_dp-y**2+(1.0_dp/3.0_dp)*y**2*log(y)
        t = v + (4.0/3.0)*y**2
        
        j_fn = a*f**2/w * exp(-v*b*w**(1.5_dp) / f)/t**2
        
    end function j_fn

    function j_q(ts, f, ef, w)
        real(dp) :: ts, f, ef, w
        real(dp) :: j_q
        
        integer, parameter :: n = 300
        integer :: i
        real(dp) :: h, y(n), e
        
        h = 30.0_dp/(n-1)
        do i = 1, n
            e = (i-1)*h
            y(i) = fermi(ef, w, ts, e)*d(ef, w, f, e, 1000)*h
        end do
        
        j_q = c_zs*sum(y)
        
    end function j_q

    function fermi(ef, w, ts, e)
        real(dp) :: ef, w, ts, e
        real(dp) :: fermi

        real(dp) :: a, b
        
        a = c_kb*ts
        b = (e-ef)/a
        
        fermi = a*log(1.0 + exp(-b))
        
    end function fermi

    function d(ef, w, f, e, n)
        real(dp) :: ef, w, f, e
        real(dp) :: d
        integer :: i, n
        
        real(dp) :: u(n)
        real(dp) :: xc, xmax, alpha, k, t1, r1
        real(dp) :: h, x
		complex(16) :: r, phi, j = (0,1.0)
        
        call find_xc_xmax(ef, w, f, xc, xmax)
		
        alpha = 2.0d-18/c_hp*c_me/c_hp*c_qe

        k = sqrt(alpha*e)

        h = (xmax-xc)/(n-1.0)

        do i = 1, n
            x = xc + (i-1)*h
            u(i) = 0.5*alpha*h**2*(e-v(ef, w, f, x, xc, xmax))-1.0
        end do

        r = -1.0 / (u(n) + j*k*h)

        do i = n-1, 2, -1
            r = -1.0 / (2*u(i) + r)
        end do

        phi = 2*j*k*h / (r + j*k*h + u(1))

        d = 1-abs(phi-1.0)**2
    end function d

    subroutine find_xc_xmax(ef, w, f, xc, xmax)
        real(dp), intent(in)  :: ef, w, f
        real(dp), intent(out) :: xc, xmax
        
        real(dp) :: a, b, c, dd

        c = ef+w
        a = c_cs**2/(4*c)
        b = f/c

        dd = 1.0_dp-4*a*b

        xc = (1.0_dp-sqrt(dd))/(2*b)
        xmax = (1.0_dp+sqrt(dd))/(2*b)

    end subroutine find_xc_xmax

    function v(ef, w, f, x, a, b)
        real(dp) :: f, ef, w, x, v, a, b

        if ((x<=a) .or. (x>=b)) then
            v = 0.0_dp
        else
            v = ef + w - c_cs**2/(4*x) - f*x
        endif

    end function v

    function v_pit(e0, x, a, b)
        real(dp) :: v_pit, e0, x, a, b

        if ((x<=a) .or. (x>=b)) then
            v_pit = 0.0_dp
        else
            v_pit = e0
        endif

    end function v_pit

    function d_pit(e0, e, n)
        real(dp) :: d_pit, e0, e
        integer :: i, n
        complex(16) :: j = (0,1.0)
        real(dp) :: u(n)
        complex(16) :: r, phi
        real(dp) :: xc, xmax, alpha, k, t1, r1
        real(dp) :: h, x, l

        xc = 0.0_dp
        xmax = l

        alpha = 2.0d-18/c_hp*c_me/c_hp*c_qe

        k = sqrt(alpha*e)

        h = (xmax-xc)/(n-1.0)

        do i = 1, n
            x = xc + (i-1)*h
            u(i) = 0.5*alpha*h**2*(e-v_pit(e0, x, xc, xmax))-1.0
        end do

        r = -1.0 / (u(n) + j*k*h)

        do i = n-1, 2, -1
            r = -1.0 / (2*u(i) + r)
        end do

        phi = 2*j*k*h / (r + j*k*h + u(1))

        d_pit = 1-abs(phi-1.0)**2
    end function d_pit

end module emission
program shrod
    use emission
	implicit none



!    open (1, file='solution.dat')
!    do while (F<14)
!        j1 = J(Ts, F)
!        j2 = J_FN(Ts, F)+J_RDS(Ts, F)
!        write (*, *), F, abs(j1-j2)/abs(j2)
!        F = F + 0.1
!    end do

!    double precision :: E, D_exact, E0, D_calc
!    double precision :: L = 10.0, alpha
!    open (1, file='solution.dat')
!    E = 0.1
!    E0 = -10.0
!    alpha = 2*me*qe*1d-18/hp**2
!    do while (E<10)
!        D_exact = 4*E*(E+abs(E0))/(4*E*(E+abs(E0))+E0**2*sin(L*sqrt(alpha*(E+abs(E0))))**2)
!        D_calc = D_pit(E0, E, 500)
!        write (1, *), E, D_calc, D_exact
!        print *, abs(D_calc - D_exact)/D_exact
!        E = E+0.1
!    end do

end program shrod

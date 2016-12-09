    module global_data

      implicit none

      integer(4),       public, parameter   :: rk = selected_real_kind(8)
      integer(4),       public, parameter   :: npolar = 16
      integer(4),       public, parameter   :: nalpha = 200
      integer(4),       public, parameter   :: pmax   = 6
      integer(4),       public, parameter   :: cmax   = 4
      integer(4),       public, parameter   :: smax   = 200
      real(kind=rk),    public              :: pi
      real(kind=rk),    public              :: mu(npolar)
      real(kind=rk),    public              :: wgt(npolar)
      real(kind=rk),    public              :: c(1:cmax)=(/0.8_rk,0.9_rk,0.99_rk,0.999999_rk/)
      real(kind=rk),    public              :: rho(smax,cmax,pmax)
      complex(kind=rk), public, allocatable :: um(:,:)

    end module global_data

    program fourier_analysis

      use global_data, only: rk, pi

      implicit none

      integer(4),    parameter :: sweeps=1

      pi=2.0_rk*acos(0.0_rk)

      call spectral_radius(sweeps)
      call plot

    contains

    subroutine spectral_radius(sweeps)

      use global_data, only: rk, pmax, cmax, smax, c, rho

      implicit none

      integer(4)                 :: p
      integer(4)                 :: r
      integer(4)                 :: s
      integer(4),    intent(in)  :: sweeps
      real(kind=rk)              :: sigt_h
      real(kind=rk)              :: del(smax,pmax)

      call set_quadrature
      call set_del(del)

      rho=0.0_rk
      do p=1,pmax
        do r=1,cmax
          do s=1,smax
            sigt_h=del(s,p)/real(p,rk)
            call get_rho(p, sweeps, del(s,p), c(r), sigt_h, rho(s,r,p))
          enddo
        enddo
      enddo

    end subroutine spectral_radius

    subroutine set_quadrature

      use global_data, only: rk, npolar, pi, mu, wgt

      implicit none

      integer(4)                :: i
      integer(4)                :: j
      integer(4)                :: m
      real(kind=rk)             :: p1
      real(kind=rk)             :: p2
      real(kind=rk)             :: p3
      real(kind=rk)             :: pp
      real(kind=rk)             :: xm
      real(kind=rk)             :: xl
      real(kind=rk)             :: z
      real(kind=rk)             :: z1
      real(kind=rk)             :: eps

      mu=0.0_rk
      wgt=0.0_rk
      eps=epsilon(0.0_rk)

      m=(npolar+1)/2
      xm=0.0_rk
      xl=1.0_rk

      do i=1,m
        z=cos(pi*(real(i,rk)-0.25_rk)/(real(npolar,rk)+0.5_rk))
    1   continue
        p1=1.0_rk
        p2=0.0_rk
        do j=1,npolar
          p3=p2
          p2=p1
          p1=((2.0_rk*real(j,rk)-1.0_rk)*z*p2-(real(j,rk)-1.0_rk)*p3)/real(j,rk)
        enddo
        pp=real(npolar,rk)*(z*p1-p2)/(z*z-1.0_rk)
        z1=z
        z=z1-p1/pp
        if (abs(z-z1) > eps) go to 1
        mu(npolar+1-i) =xm-xl*z
        mu(i)          =xm+xl*z
        wgt(npolar+1-i)=2.0_rk*xl/((1.0_rk-z*z)*pp*pp)
        wgt(i)         =wgt(npolar+1-i)
      enddo

    end subroutine set_quadrature

    subroutine set_del(del)

      use global_data, only: rk, pmax, smax

      implicit none

      integer(4)                 :: p
      integer(4)                 :: s
      integer(4),    parameter   :: sm1=100
      real(kind=rk)              :: ds
      real(kind=rk), intent(out) :: del(smax,pmax)

      del=0.0_rk
      do p=1,pmax
        ds=log(100.0_rk)/real(smax-sm1,rk)
        do s=1,smax 
          if (s <= sm1) then
            del(s,p)=(1.0_rk/real(sm1,rk))*real(s,rk)
          else
            del(s,p)=exp(ds*real(s-sm1,rk))
          endif
        enddo
      enddo

    end subroutine set_del

    subroutine get_rho(p, sweeps, del, c, sigt_h, rho)

      use global_data, only: rk, nalpha, pi

      implicit none

      integer(4)                 :: i
      integer(4),    intent(in)  :: p
      integer(4),    intent(in)  :: sweeps
      real(kind=rk)              :: alpha
      real(kind=rk)              :: omega
      real(kind=rk), intent(in)  :: del
      real(kind=rk), intent(in)  :: sigt_h
      real(kind=rk), intent(in)  :: c
      real(kind=rk), intent(out) :: rho

      alpha=0.0_rk
      do i=1,nalpha
        alpha=(pi/del)*(real(i,rk)/real(nalpha,rk))
        call buildm(p, sweeps, sigt_h, c, del, alpha)
        call get_eigenvalue(p,omega)
        rho=max(rho,omega)
      enddo

    end subroutine get_rho

    subroutine buildm(p, sweeps, sigt_h, c, del, alpha)

      use global_data, only: rk, npolar, mu, um

      implicit none

      integer(4)                   :: n
      integer(4),       intent(in) :: p
      integer(4),       intent(in) :: sweeps
      real(kind=rk)                :: d
      real(kind=rk)                :: s
      real(kind=rk)                :: diffco
      real(kind=rk)                :: fhat
      real(kind=rk)                :: beta
      real(kind=rk)                :: tau
      real(kind=rk),    intent(in) :: sigt_h
      real(kind=rk),    intent(in) :: c
      real(kind=rk),    intent(in) :: del
      real(kind=rk),    intent(in) :: alpha
      complex(kind=rk)             :: am(p,p)
      complex(kind=rk)             :: bm(p,p)
      complex(kind=rk)             :: ab(p,p)

      if (allocated(um)) deallocate(um)
      allocate(um(p,p))
      um=(0.0_rk,0.0_rk)

      do n=1,npolar/2
        tau=sigt_h/mu(n)
        call get_beta(tau, beta)
        d=(1.0_rk-beta)/2.0_rk
        s=(1.0_rk+beta)/2.0_rk
        call setmat(p, alpha, del, d, s, am)
        d=-1.0_rk/tau
        s= 1.0_rk/tau
        call setmat(p, alpha, del, d, s, bm)
        bm=am+bm
        call invcmat(p, bm)
        ab=MATMUL(am,bm)
        call addum(n, p, c, ab)
      enddo

      diffco=1.0_rk/(3.0_rk*del)
      fhat=c*sigt_h/(2.0_rk*diffco*(1.0_rk-cos(alpha*del))+(1.0_rk-c)*del)

      call fnlum(p, sweeps, fhat)

    end subroutine buildm

    subroutine get_beta(tau, beta)

      use global_data, only: rk

      implicit none

      real(kind=rk)              :: tau3
      real(kind=rk)              :: tau5
      real(kind=rk)              :: tau7
      real(kind=rk), intent(in)  :: tau
      real(kind=rk), intent(out) :: beta

      if (tau < 0.01_rk) then
        tau3=tau *tau*tau
        tau5=tau3*tau*tau
        tau7=tau5*tau*tau
        beta=tau/6.0_rk-tau3/360.0_rk+tau5/15120.0_rk-tau7/604800.0_rk
      else
        beta=1.0_rk/tanh(tau/2.0_rk)-2.0_rk/tau
      endif

    end subroutine get_beta

    subroutine setmat(m, alpha, del, d, s, a)

      use global_data, only: rk

      implicit none

      integer(4)                    :: i
      integer(4),       intent(in)  :: m
      real(kind=rk)                 :: re
      real(kind=rk)                 :: im
      real(kind=rk),    intent(in)  :: alpha
      real(kind=rk),    intent(in)  :: del
      real(kind=rk),    intent(in)  :: d
      real(kind=rk),    intent(in)  :: s
      complex(kind=rk), intent(out) :: a(m,m)

      a=(0.0_rk,0.0_rk)
      do i=1,m
        a(i,i)=d
      enddo
      do i=1,m-1
        a(i,i+1)=s
      enddo
      re=s*cos(alpha*del)
      im=s*sin(alpha*del)
      a(m,1)=a(m,1)+cmplx(re,im,rk)

    end subroutine setmat

    subroutine invcmat(m, a)

      use global_data, only: rk

      implicit none

      integer(4)                      :: info
      integer(4),       intent(in)    :: m
      integer(4)                      :: ipiv(m)
      complex(kind=rk), intent(inout) :: a(m,m)
      complex(kind=rk)                :: w(m)

      ipiv=0
      call zgetrf(m, m, a, m, ipiv, info)
      if (info /= 0) stop ' Trouble with LU factorization.'
      w=(0.0_rk,0.0_rk)
      call zgetri(m, a, m, ipiv, w, m, info)
      if (info /= 0) stop ' Trouble with finding matrix inverse.'

    end subroutine invcmat

    subroutine addum(n, m, c, ab)

      use global_data, only: rk, wgt, um

      implicit none

      integer(4)                   :: i
      integer(4)                   :: j
      integer(4),       intent(in) :: n
      integer(4),       intent(in) :: m
      real(kind=rk),    intent(in) :: c
      complex(kind=rk), intent(in) :: ab(m,m)

      do i=1,m
        do j=1,m
          um(i,j)=um(i,j)+c*wgt(n)*ab(i,j)
        enddo
      enddo

    end subroutine addum

    subroutine fnlum(m, sweeps, fhat)

      use global_data, only: rk, um

      implicit none

      integer(4)                   :: i,j
      integer(4)                   :: s
      integer(4),       intent(in) :: m
      integer(4),       intent(in) :: sweeps
      real(kind=rk),    intent(in) :: fhat
      complex(kind=rk)             :: um1(m,m)
      complex(kind=rk)             :: ums(m,m)
      complex(kind=rk)             :: fm (m,m)

      um1=(1.0_rk,0.0_rk)
      do i=1,m
        do j=1,m
          um(i,j)=um(i,j)+fhat*m*(um(i,j)-1.0_rk)
        enddo
      enddo
      !um=um+fhat*(um-1.0_rk)

      !ums=um
      !do s=1,sweeps-1
      !  if (s == sweeps-1) um1=ums
      !  ums=MATMUL(um,ums)
      !enddo
      !fm=fhat !*(1.0_rk,0.0_rk)
      !um=ums+MATMUL(fm,ums)-MATMUL(fm,um1)

    end subroutine fnlum

    subroutine get_eigenvalue(m,omega)

      use global_data, only: rk, um

      implicit none

      integer(4),       intent(in)  :: m
      integer(4)                    :: info
      real(kind=rk)                 :: rw(2*m)
      real(kind=rk),    intent(out) :: omega
      complex(kind=rk)              :: vl
      complex(kind=rk)              :: vr
      complex(kind=rk)              :: e(m)
      complex(kind=rk)              :: w(2*m)
      character(1),     parameter   :: jobvl='N'
      character(1),     parameter   :: jobvr='N'

      call zgeev(jobvl, jobvr, m, um, m, e, vl, 1, vr, 1, w, 2*m, rw, info)
      if (info /= 0) stop ' Trouble with finding matrix eigenvalues.'
      omega=maxval(abs(real(e)))

    end subroutine get_eigenvalue

    subroutine plot

      use global_data, only: rk, pmax, cmax, smax, rho

      implicit none

      integer(4)                :: i
      integer(4)                :: j
      integer(4)                :: k
      integer(4)                :: dun=1
      integer(4)                :: pun=2
      real(kind=rk)             :: del(smax,pmax)
      character(1)              :: p
      character(2)              :: cols
      character(13)             :: datafile
      character(9)              :: plotfile
      character(10)             :: postfile
      character(132)            :: datfmt 

      call set_del(del)

      write(cols,'(i2)') 1+cmax
      write(datfmt,'(a)') '(' // trim(adjustl(cols)) // '(es12.5))'

      do k=1,pmax
        write(p,'(i1)') k
        write(datafile,'(a)') 'result-p' // trim(adjustl(p)) // '.dat'
        open(unit=dun,file=datafile,action='write',status='unknown')
        do i=1,smax
          write(dun,datfmt) del(i,k),(rho(i,j,k),j=1,cmax)
        enddo
        close(dun)
        write(plotfile,'(a)') 'plot-p' // trim(adjustl(p)) // '.p'
        write(postfile,'(a)') 'plot-p' // trim(adjustl(p)) // '.ps'
        open(unit=pun,file=plotfile,action='write',status='unknown')
        write(pun,'(a)') '# Gnuplot script file for plotting data'
        write(pun,'(a)') 'set autoscale                          # scale axes automatically'
        write(pun,'(a)') 'unset logscale; set logscale x         # set x log-scaling'
        write(pun,'(a)') 'unset label                            # remove any previous labels'
        write(pun,'(a)') 'set xtic auto                          # set xtics automatically'
        write(pun,'(a)') 'set ytic auto                          # set ytics automatically'
        write(pun,'(a)') 'set title "Step Characteristics, p = 1"'
        write(pun,'(a)') 'set xlabel "{/Symbol D} (mfp)" enhanced'
        write(pun,'(a)') 'set ylabel "{/Symbol r}" enhanced'
        write(pun,'(a)') '# set key 0.01,100'
        write(pun,'(a)') 'set yr [0:1]'
        write(pun,'(a)') 'plot    "' // datafile // '" using 1:2 title "c=0.8"  with linespoints , \'
        write(pun,'(a)') '        "' // datafile // '" using 1:3 title "c=0.9"  with linespoints , \'
        write(pun,'(a)') '        "' // datafile // '" using 1:4 title "c=0.99" with linespoints , \'
        write(pun,'(a)') '        "' // datafile // '" using 1:5 title "c=1.00" with linespoints'
        write(pun,'(a)') 'set size 1.0, 0.6'
        write(pun,'(a)') 'set terminal postscript portrait enhanced color dashed lw 1 "Helvetica" 14'
        write(pun,'(a)') 'set output "' // postfile  // '"'
        write(pun,'(a)') 'replot'
        write(pun,'(a)') 'set terminal x11'
        write(pun,'(a)') 'set size 1,1'
        close(pun)
      enddo

    end subroutine plot

    end program fourier_analysis

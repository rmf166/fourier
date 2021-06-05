    module global_data

      implicit none

      integer(4),       parameter   :: rk = kind(1.0d0)
      integer(4),       parameter   :: npolar = 16
      integer(4),       parameter   :: nalpha = 1000
      integer(4),       parameter   :: pmax   = 6
      integer(4),       parameter   :: cmax   = 4
      integer(4),       parameter   :: jmax   = 500
      integer(4),       parameter   :: maxs   = 4
      integer(4)                    :: sweeps
      integer(4)                    :: umax
      logical(4)                    :: linpro
      real(kind=rk)                 :: pi
      real(kind=rk)                 :: mu(npolar)
      real(kind=rk)                 :: wgt(npolar)
      real(kind=rk)                 :: c(1:cmax)=(/0.8_rk,0.9_rk,0.99_rk,0.9999_rk/)
      real(kind=rk)                 :: rho(jmax,cmax,pmax)
      character(2)                  :: sol
      complex(kind=rk), allocatable :: um(:,:)

    end module global_data

    program fourier_analysis

      use global_data, only: rk, maxs, sweeps, linpro, pi, sol

      implicit none

      integer(4)   :: i
      integer(4)   :: j
      integer(4)   :: k
      real(4)      :: start
      real(4)      :: finish
      character(2) :: solopt(1:4)=(/'DD','SC','LD','LC'/)

      pi=2.0_rk*acos(0.0_rk)

      do k=1,2
        linpro=.false.
        if (k == 2) linpro=.true.
        do j=1,4
          sol=solopt(j)
          do i=1,maxs
            sweeps=i
            call cpu_time(start)
            call spectral_radius
            call plot
            call cpu_time(finish)
            write(0,'(3a,i3,a,f6.3,a)') ' Finished Fourier Analysis for ', sol,   &
                                        ' discretization and ', i, ' sweeps in ', &
                                          (finish-start)/60.0,' minutes.'
          enddo
        enddo
      enddo

    end program fourier_analysis

    subroutine spectral_radius

      use global_data, only: rk, pmax, cmax, jmax, c, rho

      implicit none

      integer(4)                 :: p
      integer(4)                 :: r
      integer(4)                 :: s
      real(kind=rk)              :: sigt_h(jmax,pmax)

      call set_quadrature
      call set_sigt(sigt_h)

      rho=0.0_rk
      do p=1,pmax
        do r=1,cmax
          do s=1,jmax
            call get_rho(p, c(r), sigt_h(s,p), rho(s,r,p))
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

    subroutine set_sigt(sigt_h)

      use global_data, only: rk, pmax, jmax

      implicit none

      integer(4)                 :: p
      integer(4)                 :: s
      integer(4),    parameter   :: sm1=100
      real(kind=rk)              :: ds
      real(kind=rk), intent(out) :: sigt_h(jmax,pmax)

      if (sm1 >= jmax) stop ' Trouble in set_sigt, sm1 >= jmax.'
      sigt_h=0.0_rk
      do p=1,pmax
        ds=log(100.0_rk)/real(jmax-sm1,rk)
        do s=1,jmax 
          if (s <= sm1) then
            sigt_h(s,p)=(1.0_rk/real(sm1,rk))*real(s,rk)
          else
            sigt_h(s,p)=exp(ds*real(s-sm1,rk))
          endif
        enddo
      enddo

    end subroutine set_sigt

    subroutine get_rho(p, c, sigt_h, rho)

      use global_data, only: rk, nalpha, sol, umax, pi

      implicit none

      integer(4)                 :: i
      integer(4),    intent(in)  :: p
      real(kind=rk)              :: alpha
      real(kind=rk)              :: eps
      real(kind=rk)              :: omega
      real(kind=rk), intent(in)  :: sigt_h
      real(kind=rk), intent(in)  :: c
      real(kind=rk), intent(out) :: rho

      if (sol == 'DD' .or. sol == 'SC') then
        umax=p
      elseif (sol == 'LD' .or. sol == 'LC') then
        umax=2*p
      endif

      eps=0.0000000001_rk
      alpha=eps
      do i=1,nalpha+1
        call buildm(p, umax, c, sigt_h, alpha)
        call get_eigenvalue(umax,omega)
        rho=max(rho,omega)
        alpha=i*((0.5_rk*pi-eps)/nalpha)
      enddo

    end subroutine get_rho

    subroutine buildm(p, umax, c, sigt_h, alpha)

      use global_data, only: rk, npolar, linpro, mu

      implicit none

      integer(4)                   :: i
      integer(4)                   :: j
      integer(4)                   :: n
      integer(4),       intent(in) :: p
      integer(4),       intent(in) :: umax
      real(kind=rk)                :: diffco
      real(kind=rk)                :: fhat
      real(kind=rk)                :: beta
      real(kind=rk)                :: eta
      real(kind=rk),    intent(in) :: c
      real(kind=rk),    intent(in) :: sigt_h
      real(kind=rk),    intent(in) :: alpha
      complex(kind=rk)             :: arg
      complex(kind=rk)             :: expp
      complex(kind=rk)             :: expm
      complex(kind=rk)             :: fact
      complex(kind=rk)             :: d
      complex(kind=rk)             :: s
      complex(kind=rk)             :: am(p,p)
      complex(kind=rk)             :: bm(p,p)
      complex(kind=rk)             :: cm(p,p)
      complex(kind=rk)             :: uv(umax,umax)

      am =(0.0_rk,0.0_rk)
      bm =(0.0_rk,0.0_rk)
      cm =(0.0_rk,0.0_rk)
      uv =(0.0_rk,0.0_rk)

      arg=cmplx(0.0_rk,alpha,rk)
      expp=exp(arg)
      expm=exp(-arg)

      call allocum(umax)

      do n=1,npolar
        eta=sigt_h/mu(n)
        call get_beta(eta, beta)
        d=(1.0_rk-beta)/2.0_rk*expm
        s=(1.0_rk+beta)/2.0_rk*expp
        call setmat(p, d, s, am)
        call invcmat(p, am)
        d=-1.0_rk/eta*expm
        s= 1.0_rk/eta*expp
        call setmat(p, d, s, bm)
        d= 3.0_rk/eta*expm
        s= 3.0_rk/eta*expp
        call setmat(p, d, s, cm)
        call fnlmat(n, p, umax, c, eta, beta, am, bm, cm)
      enddo

      do j=1,p
        expp=exp( arg*(2.0_rk*real(j,rk)-1.0_rk-real(p,rk)))
        do i=1,p
          fact=(1.0_rk,0.0_rk)
          if (linpro) fact=((2.0_rk*real(p,rk)-2.0_rk*real(i,rk)+1.0_rk)*exp(-arg*2.0_rk*real(p,rk)) +&
                            (2.0_rk*real(i,rk)                  -1.0_rk)*exp( arg*2.0_rk*real(p,rk)))/&
                            (4.0_rk*real(p,rk))+0.5_rk
          expm=exp(-arg*(2.0_rk*real(i,rk)-1.0_rk-real(p,rk)))*fact
          uv(i,j)=expm*expp
        enddo
      enddo

      diffco=1.0_rk/(3.0_rk*sigt_h*real(p,rk))
      fhat=c*sigt_h/(2.0_rk*diffco*(1.0_rk-cos(2.0_rk*alpha*real(p,rk)))+(1.0_rk-c)*sigt_h*real(p,rk))

      call fnlum(umax, uv, fhat)

    end subroutine buildm

    subroutine allocum(umax)

      use global_data, only: rk, um

      implicit none

      integer(4),    intent(in) :: umax

      if (allocated(um)) deallocate(um)
      allocate(um(umax,umax))
      um=(0.0_rk,0.0_rk)

    end subroutine allocum

    subroutine get_beta(eta, beta)

      use global_data, only: rk, sol

      implicit none

      real(kind=rk), intent(in)  :: eta
      real(kind=rk), intent(out) :: beta

      if (sol == 'DD') then
        beta=0.0_rk
      elseif (sol == 'LD') then
        beta=eta/abs(eta)
      else
        beta=1.0_rk/tanh(eta/2.0_rk)-2.0_rk/eta
        if (sol == 'LC') beta=eta*beta/(eta-6.0_rk*beta)
      endif

    end subroutine get_beta

    subroutine fnlmat(n, p, umax, c, eta, beta, am, bm, cm)

      use global_data, only: rk, sol

      implicit none
      integer(4)                   :: i
      integer(4),       intent(in) :: n
      integer(4),       intent(in) :: p
      integer(4),       intent(in) :: umax
      real(kind=rk),    intent(in) :: c
      real(kind=rk),    intent(in) :: eta
      real(kind=rk),    intent(in) :: beta
      complex(kind=rk)             :: eye(p,p)
      complex(kind=rk)             :: mm(umax,umax)
      complex(kind=rk), intent(in) :: am(p,p)
      complex(kind=rk), intent(in) :: bm(p,p)
      complex(kind=rk), intent(in) :: cm(p,p)

      eye=(0.0_rk,0.0_rk)
      do i=1,p
        eye(i,i)=(1.0_rk,0.0_rk)
      enddo 

      if (sol == 'DD' .or. sol == 'SC') then
        mm=MATMUL(bm,am)
        mm=mm+eye
      elseif (sol == 'LD' .or. sol == 'LC') then
        mm(    1:p,    1:p)=MATMUL(bm,am)+eye
        mm(    1:p,p+1:2*p)=beta*MATMUL(bm,am)
        mm(p+1:2*p,    1:p)=MATMUL(cm,am)-6.0_rk/eta*eye
        mm(p+1:2*p,p+1:2*p)=beta*MATMUL(cm,am)+eye
      endif
      call invcmat(umax, mm)
      call addum(n, umax, c, mm)

    end subroutine fnlmat

    subroutine setmat(m, d, s, a)

      use global_data, only: rk

      implicit none

      integer(4)                    :: i
      integer(4),       intent(in)  :: m
      complex(kind=rk), intent(in)  :: d
      complex(kind=rk), intent(in)  :: s
      complex(kind=rk), intent(out) :: a(m,m)

      a=(0.0_rk,0.0_rk)
      do i=1,m
        a(i,i)=d
      enddo
      do i=1,m-1
        a(i,i+1)=s
      enddo
      a(m,1)=a(m,1)+s

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

    subroutine addum(n, m, c, a)

      use global_data, only: rk, wgt, um

      implicit none

      integer(4)                   :: i
      integer(4)                   :: j
      integer(4),       intent(in) :: n
      integer(4),       intent(in) :: m
      real(kind=rk),    intent(in) :: c
      complex(kind=rk), intent(in) :: a(m,m)

      do i=1,m
        do j=1,m
          um(i,j)=um(i,j)+0.5_rk*c*wgt(n)*a(i,j)
        enddo
      enddo

    end subroutine addum

    subroutine fnlum(m, uv, fhat)

      use global_data, only: rk, sweeps, um

      implicit none

      integer(4)                   :: i
      integer(4),       intent(in) :: m
      real(kind=rk),    intent(in) :: fhat
      complex(kind=rk)             :: eye(m,m)
      complex(kind=rk)             :: um0(m,m)
      complex(kind=rk)             :: um1(m,m)
      complex(kind=rk), intent(in) :: uv (m,m)

      if (sweeps == 1) then
        eye=(0.0_rk,0.0_rk)
        do i=1,m
          eye(i,i)=(1.0_rk,0.0_rk)
        enddo
        um=um+fhat*MATMUL(uv,um-eye)
      else
        i=1
        um0=um
        do while (i < sweeps)
          um1=um
          um=MATMUL(um0,um)
          i=i+1
        enddo
        um=um+fhat*MATMUL(uv,um-um1)
      endif

    end subroutine fnlum

    subroutine get_eigenvalue(m,omega)

      use global_data, only: rk, um

      implicit none

      integer(4)                    :: i
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
      omega=0.0_rk
      do i=1,m
        omega=max(omega,abs(e(i)))
      enddo

    end subroutine get_eigenvalue

    subroutine plot

      use global_data, only: rk, pmax, cmax, jmax, sweeps, linpro, rho, sol

      implicit none

      integer(4)                :: i
      integer(4)                :: j
      integer(4)                :: k
      integer(4)                :: dun=1
      real(kind=rk)             :: sigt_h(jmax,pmax)
      character(1)              :: p
      character(1)              :: s
      character(2)              :: cols
      character(10)             :: caserun
      character(20)             :: datafile
      character(132)            :: datfmt 

      call set_sigt(sigt_h)

      write(cols,'(i2)') 1+cmax
      write(datfmt,'(a)') '(' // trim(adjustl(cols)) // '(es12.5))'

      do k=1,pmax
        write(p,'(i1)') k
        write(s,'(i1)') sweeps
        if (linpro) then
          write(caserun,'(a)') '-lp' // trim(adjustl(p)) // '-s' // trim(adjustl(s)) // '-' // sol
        else
          write(caserun,'(a)') '-p' // trim(adjustl(p)) // '-s' // trim(adjustl(s)) // '-' // sol
        endif
        write(datafile,'(a)') 'result' // trim(adjustl(caserun)) // '.dat'
        open(unit=dun,file=datafile,action='write',status='unknown')
        do i=1,jmax
          write(dun,datfmt) sigt_h(i,k),(rho(i,j,k),j=1,cmax)
        enddo
        close(dun)
      enddo

    end subroutine plot

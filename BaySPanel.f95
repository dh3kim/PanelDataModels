! Bayesian estimation of (balanced) spatial panel data

Module Precision
 Integer, Parameter :: sp = Kind(0.0E0)
 Integer, Parameter :: dp = Kind(0.0D0)
 Integer, Parameter :: wp = dp
End Module Precision

MODULE DataPaths
IMPLICIT none
CHARACTER(100)::                 datapath = ".././"
CHARACTER(200)::                 input
CHARACTER(200)::                 Wmatrix_input
END MODULE DataPaths

MODULE DataParameters
USE Precision
! USE LA_PRECISION, ONLY:          WP => DP
USE F95_LAPACK, ONLY:            LA_POSV, LA_SYSV, LA_POTRF
USE F95_LAPACK, ONLY:                    LA_SYEV
IMPLICIT none
 INTEGER, PARAMETER ::            nobs = 400  
 INTEGER, PARAMETER ::            N = 100
 INTEGER, PARAMETER ::            T = 4
 INTEGER, PARAMETER ::            KVars = 2   
 INTEGER, DIMENSION(nobs) ::      CityID 
 REAL(wp), DIMENSION(nobs) ::         Y
 REAL(wp), DIMENSION(nobs,KVars) ::   X

 REAL(wp), DIMENSION(N, N) ::         W

! Used in OLS
 REAL(wp), DIMENSION(KVars)         ::  XTY
 REAL(wp), DIMENSION(KVars, KVars)  ::  XTX
 REAL(wp), DIMENSION(KVars)             ::  beta
 REAL(wp), DIMENSION(nobs)              ::  e       ! OLS residuals
 REAL(wp)                               ::  sigma_sq

 REAL(wp), DIMENSION(nobs) :: WY
 REAL(wp), DIMENSION(nobs, KVars):: WX

!Bayesian
 INTEGER, PARAMETER ::                   ndraw = 500000  ! number of draw
 INTEGER, PARAMETER ::                   nburn = 450000  ! number of burn-in

 REAL(wp), DIMENSION(KVars)             ::   beta
 REAL(wp)                               ::   sigma_sq
 REAL(wp)                               ::   rho
 REAL(wp), DIMENSION(N)                 ::   MU
 REAL(wp)                               :: sigma_mu_sq


REAL(wp), DIMENSION(KVars) ::         beta0
 REAL(wp), DIMENSION(KVars, KVars) ::  M0
 REAL(wp) ::                           s0, v0
 REAL(wp) ::                           h0, p0
 REAL(wp) ::                           rho_min, rho_max
 REAL(wp) ::                           cc
 INTEGER ::                            acc

REAL(wp), DIMENSION(N) ::              eigs 

!storage of samples
 INTEGER, PARAMETER::                     npar = KVars + 3 + N 
 REAL(wp), DIMENSION(ndraw-nburn, npar):: AllDraws

 REAL(wp), DIMENSION(nobs) ::          YTilde, ZTilde
 REAL(wp), DIMENSION(nobs,KVars) ::    XTilde

REAL(wp), DIMENSION(KVars, KVars) :: M1
REAL(wp), DIMENSION(KVars, KVars) ::  identity_M1, inverseM1
REAL(wp), DIMENSION(KVars) :: Term1
REAL(wp), DIMENSION(KVars) :: str
REAL(wp), DIMENSION(KVars) :: beta1
CHARACTER(5) :: UPLO

REAL(wp) :: p1, h1
REAL(wp):: v1, s1

REAL(WP), DIMENSION(nobs):: BQ

!INTEGER, ALLOCATABLE, DIMENSION(:):: iwt 
REAL(wp):: inversePsi, DTBQ, MU1

! Sampling rho
 REAL(wp) ::                    rhostar
 !REAL(wp), DIMENSION(ndraw) ::  rhostar_save
 INTEGER ::                 accept
 REAL(wp) ::                    lnp, lnpstar
 REAL(wp) ::                    ratio, u, lnratio
 REAL(wp) ::                    rnd
 REAL(wp), DIMENSION(ndraw) ::  acc_rate
 REAL(wp)  ::                    prob

 REAL(wp)::                     rhotemp
 REAL(wp) ::                    lnconditionalrho

 REAL(wp), DIMENSION(npar):: PostMeans, PostVars
 REAL(wp), DIMENSION(npar):: NSE, RNE, CD
 REAL(wp):: frac1, frac2





END MODULE DataParameters

MODULE Random
CONTAINS

FUNCTION rnorm() RESULT( fn_val )

!   Generate a random normal deviate using the polar method.
!   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
!              normal variables', Siam Rev., vol.6, 260-264, 1964.

IMPLICIT NONE
Integer, Parameter :: sp = Kind(0.0E0)
Integer, Parameter :: dp = Kind(0.0D0)
Integer, Parameter :: wp = dp

REAL(wp)  :: fn_val

! Local variables

REAL(wp)            :: u, sum
REAL(wp), SAVE      :: v, sln
LOGICAL, SAVE   :: second = .FALSE.
REAL(wp), PARAMETER :: one = 1.0_wp, vsmall = TINY( one )

IF (second) THEN
! If second, use the second random number generated on last call

  second = .false.
  fn_val = v*sln

ELSE
! First call; generate a pair of random normals

  second = .true.
  DO
    CALL RANDOM_NUMBER( u )
    CALL RANDOM_NUMBER( v )
    u = SCALE( u, 1 ) - one
    v = SCALE( v, 1 ) - one
    sum = u*u + v*v + vsmall         ! vsmall added to prevent LOG(zero) / zero
    IF(sum < one) EXIT
  END DO
  sln = DSQRT(- SCALE( LOG(sum), 1 ) / sum)
  fn_val = u*sln
END IF

RETURN
END FUNCTION rnorm

FUNCTION random_gamma1(s, first) RESULT(fn_val)
!************************************************************
! I change random Normal generator
! x = rnorm() instead of x = random_normal()
!***********************************************************

! Uses the algorithm in
! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
! gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.

! Generates a random gamma deviate for shape parameter s >= 1.

Integer, Parameter :: sp = Kind(0.0E0)
Integer, Parameter :: dp = Kind(0.0D0)
Integer, Parameter :: wp = dp

REAL(wp), INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first
REAL(wp)                :: fn_val

! Local variables
REAL(wp), SAVE  :: c, d
REAL(wp)        :: u, v, x

IF (first) THEN
  d = s -1.0_wp/3.0_wp
  c = 1.0_wp/DSQRT(9.0_wp*d)
END IF

! Start of main loop
DO

! Generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.

  DO
    x = rnorm()
    v = (1.0_wp + c*x)**3.0
    IF (v > 0.0_wp) EXIT

  END DO

! Generate uniform variable U

  CALL RANDOM_NUMBER(u)
  IF (u < 1.0_wp - 0.0331_wp*x**4) THEN
    fn_val = d*v
    EXIT
  ELSE IF (LOG(u) < half*x**2 + d*(1.0_wp - v + LOG(v))) THEN
    fn_val = d*v
    EXIT
  END IF
END DO

RETURN
END FUNCTION random_gamma1

END MODULE Random


PROGRAM BaySPanel
USE Precision 
USE DataPaths
USE DataParameters
USE Random
IMPLICIT none
INTEGER ::                                 i, ios, j
INTEGER ::                                 iter
LOGICAL ::                                 First
REAL(wp), DIMENSION(T)::                      temp
REAL(wp), DIMENSION(N)::                   Meanetilde
!*****************************************************************
!*****************************************************************
!*****************************************************************


 input = adjustl(trim(datapath))//"VAR3.DAT"
 Wmatrix_input = adjustl(trim(datapath))//"mmat.raw"

! import data
 OPEN (Unit = 8, File = input, ACTION="READ")
 DO i = 1, nobs
   READ (Unit = 8, Fmt = *, IOSTAT=ios)  Y(i), X(i,:) 
 END DO
 CLOSE (Unit = 8)
 print *, "IO status for Y, X data is ",ios

! import spatial Weights matrix
 OPEN (Unit = 8, File = Wmatrix_input, ACTION="READ")
 DO i = 1, N
   READ (Unit = 8, Fmt = *, IOSTAT=ios) W(i, :) 
 END DO
 CLOSE (Unit = 8)
 print *, "IO status for weight matrix is ",ios

!***************************************************
! Get OLS  results                                 *
!***************************************************
 CALL OLS

 CALL SETUP_PANELGM

! create numeric id for each cross-section unit
 i = 1
 do while (i <= nobs)
    j = 0
    do while( j < 100 .and. i <= nobs)
       j = j+1
       CityID(i) = j
       i = i+1
    end do
  end do


!***************************************************************
! HyperParameters                                              *
!***************************************************************
 CALL HyperParameters

print *, " "
print *, "   ************************************************   "
print *, "   * The number of draws is ", ndraw
print *, "   * The average runtime is 39 sec. per 100 draws"
print *, "   ************************************************   "
print *, "Sampling starts ..."

! inital values
  CALL  OLS 
!  sigma_sq = 100.0_wp
  rho = 0.1_wp
  sigma_mu_sq = 100.0_wp
!  do i = 1, N
!   MU(i) = rnorm()*dsqrt(sigma_mu_sq)
!  end do
MU = 10.0_wp

identity_M1 = 0.0_wp
do i = 1, KVars
 identity_M1(i,i) = 1.0_wp
enddo

First=.true.

iter = 1
CALL RANDOM_SEED
DO WHILE (iter <= ndraw)
if (iter == 100) write(*,*) iter
if (iter == 1000) write(*,*) iter
if (iter == 2000) write(*,*) iter
if (iter == 5000) write(*,*) iter
 

 !****************************************
 ! sampling beta
 !****************************************
 XTilde = X - rho*WX
 YTilde = Y - rho*WY

! Here, YTilde = ZTilde
 DO i = 1, T 
  ZTilde( (i-1)*N+1: i*N ) = YTilde( (i-1)*N+1:i*N ) - MU
 END DO

DO i=1, KVars
  DO j=1, KVars
     XTX(i,j)=dot_product(XTilde(:,i), XTilde(:,j))
  END DO 
  XTY(i)=dot_product(XTilde(:,i), ZTilde)
END DO

  M1 = M0 + XTX/sigma_sq 
  Term1 = matmul(M0, beta0) + XTY/sigma_sq

 !invert M1
  inverseM1 = identity_M1
  call la_posv(M1, inverseM1)

  beta1 = matmul(inverseM1, Term1)


! generate random sample from N(0,1)
  do i = 1, KVars
     str(i) = rnorm()
  end do

 ! cholesky decomposiiton of inverseM1 using la_potrf
  ! make upper part zero
    call la_potrf(inverseM1, UPLO = 'L')
    do i = 1, KVars 
      do j = i+1, KVars 
        inverseM1(i,j) = 0.0_wp
      end do
    end do

    beta =  dsqrt(sigma_sq)*matmul(inverseM1, str) +  beta1
  !  write(*,*) ""
  !  write(*,*) "beta is"
  !  write(*,*) beta 

!*****************************************
 ! sampling sigma_e2
 !*****************************************
  v1 = dble(real(nobs)) + v0
  e = ZTilde - matmul(XTilde, beta)
  s1 = dot_product(e, e) + s0

 ! if (iter ==1) then 
 !        First=.true. 
 !    end if
  sigma_sq = s1/(2.0_wp * random_gamma1(v1/2.0_wp, First))
 ! First = .false.
  
! write(*,*) ""
  ! print *, "sigma_e2", sigma_sq


 !**************************************************
 ! sampling sigma_u2
 !**************************************************
  h1 = h0 + dble(real(N))
  p1 = dot_product(MU, MU) + p0
 ! print *, "p1 is ", p1
 !if (iter ==1) then 
 !        First=.true. 
 !    end if
  sigma_mu_sq = p1/(2.0_wp * random_gamma1(h1/2.0_wp, First))
 ! First = .false.
 
  !write(*,*) ""
  !print *, 'sigma_u2', sigma_mu_sq
  !print *, 'sigma_u', dsqrt(sigma_mu_sq)

!***************************************************
 ! sampling u
 !***************************************************
  BQ = YTilde - matmul(XTilde, beta)
  DO i = 1, N
   DTBQ = SUM(PACK(BQ, CityID == i))
   
   inversePsi = (sigma_mu_sq*sigma_sq)/(sigma_mu_sq*real(T) + sigma_sq)
   MU1 = (sigma_mu_sq*DTBQ)/(sigma_mu_sq*real(T) + sigma_sq)

   MU(i) = dsqrt(inversePsi)*rnorm() + MU1

  END DO
  ! write(*,*) ""
  !  print *, 'MU is',MU(1:4)

!***************************************************************** 
 ! Sampling rho with Metropolis-Hastings Algorithm                *
 !*****************************************************************
!print *, "sampling rho"
   rhostar = rho + cc * rnorm()
     accept = 0
     do while (accept == 0)
       if (( rhostar > rho_min).and.(rhostar < rho_max)) then
          accept = 1
       else
          rhostar = rho + cc * rnorm()
       end if
     end do

!   print *, "rhostar is", rhostar

    ! Calculate log of ratio
       rhotemp = rho
       call WEIGHT 
       lnp = lnconditionalrho
!       print *, "lnp is ", lnp

       rhotemp = rhostar
       call WEIGHT 
       lnpstar = lnconditionalrho
!       print *, "lnpstar is ", lnpstar

   ! accept rhostar with prob. of ratio
      call random_number(u)

     ! if ((lnpstar - lnp) > dexp(1.0_wp)) then  ! double-check the condition
       if ((lnpstar - lnp) > 0.0) then
         prob = 1.0_wp
      else
         ratio = exp(lnpstar - lnp)
         if (ratio < 1.0_wp) then
           prob = ratio
         else
           prob = 1.0_wp
         end if
     end if
!     print *, "u is ", u, "prob is ", prob
     if (u < prob) then
        rho = rhostar
        acc = acc + 1
     end if 
     
     acc_rate(iter) = real(acc)/real(iter)

    ! update cc based on std of rho draws
      if (acc_rate(iter) < 0.4_wp) then
          cc = cc/1.1_wp
       else if (acc_rate(iter) > 0.6_wp) then
          cc = cc*1.1_wp
       end if
     !print *, "acc_rate is", acc_rate(iter)
  
!print *, "rho is ", rho



   ! save samples
    if (iter > nburn) then
       AllDraws(iter-nburn, 1:KVars) = beta
       AllDraws(iter-nburn, KVars+1) = dsqrt(sigma_mu_sq)
       AllDraws(iter-nburn, KVars+2) = dsqrt(sigma_sq)
       AllDraws(iter-nburn, KVars+3) = rho
       AllDraws(iter-nburn, KVars+4:N+KVars+3) = MU 
    end if

iter = iter + 1
END DO 


CALL MC_Statistics(AllDraws, ndraw-nburn, npar, PostMeans, PostVars, NSE, RNE)
  write(*,*) "MC_Stat. done"
frac1 = 0.1_wp
frac2 = 0.5_wp

CALL Geweke_Diagnostic(AllDraws, ndraw-nburn, npar, frac1, frac2, CD)
  write(*,*) "Geweke_Diag. done"
open(unit=3, file='./result_bay_MU10', action="write")
write(3,*) "The number of draws", ndraw
write(3,*) "The number of burns", nburn
write(3,*) "acceptance rate is ", acc_rate(ndraw) 
write(3,*) "Post.Mean, Post.S.Error, NSE, RNE, CD"
do i = 1, KVars + 3 + N
write(3,*) PostMeans(i), dsqrt(PostVars(i)), NSE(i), RNE(i), cd(i)
end do



END PROGRAM BaySPanel

SUBROUTINE WEIGHT
USE Precision 
USE DataParameters
IMPLICIT none
INTEGER ::           i, j
INTEGER, DIMENSION(nobs) ::  ipivot
INTEGER, DIMENSION(KVars) :: ipiv
INTEGER ::                   infor
REAL(wp) :: logdetx, logdet
REAL(wp), DIMENSION(N) :: detp

!*****************************************************************
! Get the log of the Jacobian for this value of rho
!*****************************************************************  
!  print *, "max eigenvalues is ", maxval(eigs)
!  print *, "min eigenvalues is ", minval(eigs)
!  print *,"check1"
  detp = 1.0_wp - rhotemp*eigs
!  print *, "max detp is", maxval(detp)
!  print *, "min detp is", minval(detp)
  Where (detp <= 0.0_wp) 
     detp =0.001_wp
  end where
  logdet = real(T)*sum ( dlog(detp) )
!  print *, "logdet is ", logdet
  

!*****************************************************************
! value of log of the full conditional distribution of rho 
! under noninformative prior
! It is used to calculate acceptance probability
! double-check the use of this function 
!*****************************************************************
XTilde = X - rhotemp*WX
YTilde = Y - rhotemp*WY
DO i = 1, nobs
 ZTilde(i) = YTilde(i) - MU(CityID(i))
END DO
e = ZTilde - matmul(XTilde, beta)

!*****************************************************************
! Calculate log of determinant of XTilde'XTilde
!*****************************************************************
!  do i = 1, KVars
!   do j = 1, KVars
!     XTX(i, j) = dot_product(XTilde(:,i), XTilde(:,j))
!   end do
! end do

!  call LA_GETRF(KVars, KVars, XTX, KVars, ipiv, infor) 

!    logdetx = 0.0
!    do i = 1, KVars 
!      logdetx = logdetx + log(abs(XTX(i,i)))
!    end do

!*****************************************************************
! log of conditional distribution of rho 
!*****************************************************************
lnconditionalrho = logdet - dot_product(e, e)/(2.0_wp*sigma_sq)
             

RETURN
END SUBROUTINE WEIGHT




SUBROUTINE OLS
USE Precision
USE DataParameters
IMPLICIT none
INTEGER::          i, j
!*****************************************************************
! Solve for OLS estimator of beta, using the LU decomposition of
! X'X and right-hand side vector X'Y. Make use of the columns of
! the X matrix and dot-products for efficiency.
!*****************************************************************
  do i = 1, KVars
    do j = 1, KVars
      XTX(i,j) = dot_product(X(:,i),X(:,j))
    enddo
    XTY(i) = dot_product(X(:,i),Y)
  enddo

  call LA_POSV(XTX,XTY)
  beta = XTY
  print *, "OLS estimators for beta is", beta
 
  e = Y - matmul(X, beta)
  sigma_sq = dot_product(e, e)/real(nobs)
  print *, "OLS sigma squared is ", sigma_sq

RETURN
END SUBROUTINE OLS

SUBROUTINE SETUP_PANELGM
USE Precision
USE DataParameters
IMPLICIT none
INTEGER ::   i

  DO i = 1, T
    WX( (i-1)*N+1 : i*N, :) = MATMUL( W, X((i-1)*N+1:i*N, :) )
  END DO
  
  DO i = 1, T
    WY( (i-1)*N+1 : i*N) = MATMUL( W, Y((i-1)*N+1:i*N) )
  END DO

CALL LA_SYEV(W, eigs)
  print *, "max eigenvalues is ", maxval(eigs)
  print *, "min eigenvalues is ", minval(eigs)


RETURN
END SUBROUTINE SETUP_PANELGM


SUBROUTINE HyperParameters
USE Precision
USE DataParameters
IMPLICIT none
INTEGER ::          i

! prior
do i =1, KVars
beta0(i) = 0.0_wp
end do

M0 = 0.0_wp
do i=1, KVars
 M0(i,i) = 0.00000001_wp
end do

v0 = 5.0_wp
S0 = 500.0_wp

h0 = 5.0_wp
p0 = 500.0_wp

rho_min = 0.0_wp
rho_max = 1.0_wp/maxval(eigs)

cc = 0.2
acc = 0


RETURN
END SUBROUTINE HyperParameters


SUBROUTINE Geweke_Diagnostic(Draws, ndraws, npars, frac1, frac2, CD )
!USE Precision
USE LA_PRECISION, ONLY:    wp => dp
IMPLICIT none
INTEGER, INTENT(IN)::               ndraws, npars
REAL(wp), INTENT(IN)::              frac1, frac2
REAL(wp), DIMENSION(ndraws, npars), INTENT(IN):: Draws
REAL(wp), DIMENSION(npars), INTENT(OUT):: CD
REAL(wp), DIMENSION(npars)::        PostMeans_a, PostMeans_b
REAL(wp), DIMENSION(npars)::        PostVars_a, PostVars_b
REAL(wp), DIMENSION(npars)::        NSE_a, NSE_b
REAL(wp), DIMENSION(npars)::        RNE_a, RNE_b
INTEGER::                           ndraws_a, ndraws_b
INTEGER::                           start_a, start_b
INTEGER::                           i

!********************************************************************
! Using MCMC draws, it returns Geweke (1992)'s convergence          * 
! diagonstic, CD.                                                   *
! Arguments                                                         *
!  Draws: (ndraws * npars) matrix of MCMC draws                     *
!  ndraws: the number of draws                                      *
!  npars: the number of parameters                                  *
!  frac1: first fraction of draws used for CD                       *
!  frac2: second fraction of draws used for CD                      * 
!  CD: convergence diagonstic statistic                             *
!********************************************************************

ndraws_b = floor(frac2*dble(ndraws))
start_b =  ndraws - ndraws_b + 1
CALL MC_Statistics(Draws(start_b:ndraws,:), ndraws_b, npars, &
                   PostMeans_b, PostVars_b, NSE_b, RNE_b)

ndraws_a = floor(frac1*dble(ndraws))
start_a = 1 
CALL MC_Statistics(Draws(start_a:ndraws_a,:), ndraws_a, npars, &
                   PostMeans_a, PostVars_a, NSE_a, RNE_a)

CD = (PostMeans_a - PostMeans_b)/(NSE_a + NSE_b)

END SUBROUTINE


SUBROUTINE MC_Statistics(Draws, ndraws, npars, PostMeans, PostVars, NSE, RNE)
!USE Precision
USE LA_PRECISION, ONLY:    wp => dp
IMPLICIT none
INTEGER, INTENT(IN) ::                            ndraws, npars
REAL(wp), DIMENSION(ndraws, npars), INTENT(IN) :: Draws
REAL(wp), DIMENSION(npars), INTENT(OUT) ::        PostMeans, PostVars
REAL(wp), DIMENSION(npars), INTENT(OUT)::         NSE, RNE
REAL(wp), ALLOCATABLE, DIMENSION(:,:)::           Covars
REAL(wp), DIMENSION(npars) ::                     SG 
INTEGER ::                                        lags
INTEGER ::                                        i, j

!***************************************************************************
! Using MCMC draws, it returns posterior mean, variance, NSE (numerical    *
! Standard error, and RNE (relative numerical error).                      *
! Arguments                                                                *
!  Draws: (ndraws * npars) matrix of draws                                 *
!  ndraws: the number of draws                                             *
!  npars: the number of parameters                                         *
!  PostMeans: posterior means of the parameters                            *
!  PostVars: posterior variances of the parameters                         *
!  NSE: numerical standard error                                           *
!  RNE: relative numerical error                                           *
!***************************************************************************

! posterior means 
PostMeans = SUM(Draws,1)/dble(ndraws)

! posterior variances
do i = 1, npars
 PostVars(i) = SUM( (Draws(:,i)**2) )/dble(ndraws) - PostMeans(i)**2
end do
! variance of posterior mean is SG/ndraws
! numerical standard error (nse) is dsqrt(SG/ndraws)
! relative numerical error (rne) is PostVar/SG
lags = floor( dble(ndraws)**(0.25_wp) ) + 1
Allocate( Covars(lags, npars) )

 do i = 1, npars
  do j = 1, lags 
    Covars(j, i) = &
        dot_product( Draws(1:ndraws-j, i) - PostMeans(i), &
                     Draws(j+1:ndraws, i) - PostMeans(i) )/ &
        dble(ndraws)
  end do
 end do  

 do i = 1, npars
  SG(i) = PostVars(i) + ( 2.0_wp*dot_product( dble((/(i,i=lags,1,-1)/)), &
                       Covars(:,i)) )/ dble(lags+1) 
 end do

 RNE = PostVars/SG
 NSE = dsqrt( SG/dble(ndraws) ) 

Deallocate ( Covars )

END SUBROUTINE 



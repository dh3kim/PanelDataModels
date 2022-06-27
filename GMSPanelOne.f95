MODULE DataPaths
IMPLICIT none
CHARACTER(100)::                 saepath = "./"
CHARACTER(200)::                 input
CHARACTER(200)::                 Wmatrix_input
END MODULE DataPaths

MODULE DataParameters
 USE LA_PRECISION, ONLY:          WP => SP
 USE F95_LAPACK, ONLY:            LA_POSV, LA_SYSV, LA_POTRF
IMPLICIT none
 INTEGER, PARAMETER ::            nobs = 400  
 INTEGER, PARAMETER ::            N = 100
 INTEGER, PARAMETER ::            T = 4
 INTEGER, PARAMETER ::            KVars = 2   
 INTEGER, DIMENSION(nobs) ::      id
 REAL, DIMENSION(nobs) ::         Y
 REAL, DIMENSION(nobs,KVars) ::   X

 REAL, DIMENSION(N, N) ::         W

 ! Used in OLS
 REAL(wp), DIMENSION(KVars)         ::  XTY
 REAL(wp), DIMENSION(KVars, KVars)  ::  XTX
 REAL, DIMENSION(KVars)             ::  beta
 REAL, DIMENSION(nobs)              ::  e       ! OLS residuals
 REAL                               ::  sigma_sq

 ! Used in GM
 REAL, DIMENSION(N, N) :: I_N
 REAL, DIMENSION(T, T) :: I_T
 REAL, DIMENSION(T, T) :: J_T

 REAL, DIMENSION(nobs, nobs) :: ItW
 REAL, DIMENSION(nobs, nobs) :: Q      ! Q0 matrix
 REAL, DIMENSION(nobs, nobs) :: QQ     ! Q1 matrix

 REAL, DIMENSION(nobs) :: We, WWe
 REAL, DIMENSION(nobs) :: QQetilde

 ! moments
 REAL, DIMENSION(nobs) :: Qe, QWe, QWWe
 REAL, DIMENSION(nobs) :: QQe, QQWe, QQWWe

 REAL ::                  eQe, eQWe, eQWWe
 REAL ::                  WeQWe, WWeQWe
 REAL ::                  WeQWWe, WWeQWWe
 REAL ::                  eQQe, eQQWe, eQQWWe
 REAL ::                  WeQQWe, WWeQQWe
 REAL ::                  WeQQWWe, WWeQQWWe

 REAL ::                          trWTW

 INTEGER, PARAMETER ::            n_moments = 6
 INTEGER, PARAMETER ::            n_var = 4
 INTEGER, PARAMETER ::            n_param = 3

 INTEGER, PARAMETER ::            n_moments0 = 3
 INTEGER, PARAMETER ::            n_var0 = 3
 INTEGER, PARAMETER ::            n_param0 = 2


 REAL, DIMENSION(n_moments, n_var) ::     Gn
 REAL, DIMENSION(n_moments) ::            small_gn

 REAL, DIMENSION(n_param) ::              delta

 CHARACTER(100)  ::                       estimate

 REAL, ALLOCATABLE, DIMENSION(:) ::       delta_temp
 REAL, ALLOCATABLE, DIMENSION(:) ::       delta0
 REAL, ALLOCATABLE, DIMENSION(:) ::       Error
 REAL, ALLOCATABLE, DIMENSION(:) ::       SE

 REAL, ALLOCATABLE, DIMENSION(:,:) ::     Gn0
 REAL, ALLOCATABLE, DIMENSION(:) ::       small_gn0

 REAL, ALLOCATABLE, DIMENSION(:) ::       g
 REAL, ALLOCATABLE, DIMENSION(:,:) ::     H
 REAL, ALLOCATABLE, DIMENSION(:,:) ::     Hinverse
 REAL, ALLOCATABLE, DIMENSION(:,:) ::     identity_H
 REAL, ALLOCATABLE, DIMENSION(:) ::       ddelta

 REAL, DIMENSION(nobs) ::                 etilde
 REAL, DIMENSION(nobs) ::                 QQetilde

 REAL ::                                  sigma_one_sq

 REAL, DIMENSION(n_param0, n_param0) ::   Xitemp

 REAL, DIMENSION(N, N) ::                 WTW
 REAL, DIMENSION(N, N) ::                 WTWWTW
 REAL ::                                  trWTWWTW

 REAL, DIMENSION(N, N) ::                 WTplusW
 REAL, DIMENSION(N, N) ::                 WTWWTplusW
 REAL ::                                  trWTWWTplusW

 REAL ::                                  trWW
 REAL ::                                  trWWplusWTW

 REAL, DIMENSION(n_param, n_param) ::     Tw
 REAL, DIMENSION(n_moments, n_moments) :: Xi
 REAL, DIMENSION(n_moments, n_moments) :: identity_Xi
 REAL, DIMENSION(n_moments, n_moments) :: XiInverse
 REAL, DIMENSION(n_moments, n_var) ::     Gntemp
 REAL, DIMENSION(n_moments) ::            small_gn_temp


! Used in GLS
 REAL ::                                  rho

 REAL, DIMENSION(nobs, nobs) ::           IIrhoW
 REAL, DIMENSION(nobs) ::                 YTilde
 REAL, DIMENSION(nobs,KVars) ::           XTilde

 REAL, DIMENSION(nobs, nobs) ::           OmegInv
 REAL, DIMENSION(nobs) ::                 YYTilde
 REAL, DIMENSION(nobs,KVars) ::           XXTilde

 REAL(wp), DIMENSION(KVars) ::            XTY
 REAL(wp), DIMENSION(KVars,KVars) ::      XTX
 REAL(wp), DIMENSION(KVars,KVars) ::      identity_beta
 REAL, DIMENSION(KVars,KVars) ::          var_beta

END MODULE DataParameters



PROGRAM GMSPanel 
USE DataPaths
USE DataParameters
IMPLICIT none
INTEGER ::                                 i, ios, j
INTEGER ::                                 iter
CHARACTER :: UPLO
!*****************************************************************
! Spatial Panel Data Model                                       *
! Kappor, Kelejian, and Prucha (2004)'s GM Estimation            *
!*****************************************************************
! See line 255 to choose one between weighted GM and             *
!   partially weighted GM                                        *
!*****************************************************************

 input = adjustl(trim(saepath))//"VAR3.DAT"
 Wmatrix_input = adjustl(trim(saepath))//"mmat.raw"

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


 CALL SETUP

!***************************************************
! Get OLS  results                                 *
!***************************************************
 CALL OLS

!****************************************************
! GM estimation                                     * 
!****************************************************

 CALL SETUP_PANELGM

! call elements of sample moments                  
 CALL MOMENTS

!****************************************************
! Find inital GM estimators of rho and sigma_sq     *
!****************************************************
! initial guesses
  delta(1) = 0.1       ! initial value of rho
  delta(2) = sigma_sq

  estimate = 'initialGM' 

  allocate( delta_temp(n_param0) )
  allocate( delta0(n_param0) )
  allocate( Error(n_moments0) )
  allocate( SE(n_moments0))
  allocate( Gn0(n_moments0, n_var0) )
  allocate( small_gn0(n_moments0) )

  allocate( g(n_param0) )
  allocate( H(n_param0, n_param0) )
  allocate( Hinverse(n_param0, n_param0) )
  allocate( identity_H(n_param0, n_param0) )
  allocate( ddelta(n_param0) )

! Create an identity matrix to be used in
  identity_H = 0.0_wp
    do i = 1, n_param0
    identity_H(i,i) = 1.0_wp
  end do
 
  Gn0 = Gn(1:n_moments0, 1:n_var0)
  small_gn0 = small_gn(1:n_moments0)

  delta_temp = delta(1:n_param0)

! Newton optimization algorithm
  CALL NEWTON

  delta(1:n_param0) = delta_temp

  deallocate( delta_temp )
  deallocate( delta0 )
  deallocate( Error )
  deallocate( SE)
  deallocate( Gn0 )
  deallocate( small_gn0 )

  deallocate( g )
  deallocate( H )
  deallocate( Hinverse )
  deallocate( identity_H )
  deallocate( ddelta )

 
  print *, 'initial GM estimator of rho', delta(1)
  print *, 'initial GM estimator of sigma_sq', delta(2)

! sigma_one_sq
  etilde = e - delta(1)*We
  do i = 1, nobs
    QQetilde(i) = dot_product(QQ(i,:),etilde)
  end do
  sigma_one_sq = (1.0/real(N))*dot_product(etilde, QQetilde)
  print *, 'initial GM sigma1', sigma_one_sq


!**************************************************************
! Find the weighted GM or partially weighted estimator        *
!**************************************************************

! variance-covariance matrix of sample moments
 delta(3) = sigma_one_sq
 Xitemp = 0.0
 Xitemp(1,1) = (1.0/(real(T)-1.0))*delta(2)**2
 Xitemp(2,2) = sigma_one_sq**2

 estimate = 'PWeightedGM'  ! 'WeightedGM'

 CALL WEIGHTED_MOMENTS

 allocate( delta_temp(n_param) )
 allocate( delta0(n_param) )
 allocate( Error(n_moments) )
 allocate( SE(n_moments))
 allocate( Gn0(n_moments, n_var) )
 allocate( small_gn0(n_moments) )

 allocate( g(n_param) )
 allocate( H(n_param, n_param) )
 allocate( Hinverse(n_param, n_param) )
 allocate( identity_H(n_param, n_param) )
 allocate( ddelta(n_param) )

! Create an identity matrix to be used in
  identity_H = 0.0_wp
  do i = 1, n_param
    identity_H(i,i) = 1.0_wp
  end do

  Gn0 = Gn
  small_gn0 = small_gn
 
  delta_temp = delta

  CALL NEWTON

  delta = delta_temp

  deallocate( delta0 )
  deallocate( delta_temp )
  deallocate( Error )
  deallocate( SE)
  deallocate( Gn0 )
  deallocate( small_gn0 )

  deallocate( g )
  deallocate( H )
  deallocate( Hinverse )
  deallocate( identity_H )
  deallocate( ddelta )

print *, "GM estimator of rho", delta(1)
print *, "GM estimator of sigma_sq", delta(2)
print *, "GM estimator of sigma_one_sq", delta(3)


!***************************************************************
! Get GLS                                                      *
!***************************************************************
  rho = delta(1)
  sigma_sq = delta(2)
  sigma_one_sq = delta(3)

  CALL GLS
 
END PROGRAM GMSPanel



SUBROUTINE SETUP
USE DataParameters
IMPLICIT none
INTEGER ::          i

! Data include constant terms
 !X(:,1) = 1.0

 identity_beta = 0.0_wp
 do i = 1, KVars
   identity_beta(i,i) = 1.0_wp
 enddo

 identity_Xi = 0.0_wp
 do i = 1, n_moments
   identity_Xi(i,i) = 1.0_wp
 end do

RETURN
END SUBROUTINE SETUP



SUBROUTINE OLS
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
USE DataParameters
IMPLICIT none
INTEGER ::   i

! identity matrix of N by N
  I_N = 0.0
  do i = 1, N
    I_N(i,i) = 1.0
  end do

! identity matrix of T by T
  I_T = 0.0
  do i = 1, T
    I_T(i,i) = 1.0
  end do

! matrix of 1s
  J_T = 1.0

! kronecker product of I_T and W
  call KRON(I_T, T, T, W, N, N, ItW)

! Q0 matrix
  call KRON(I_T - J_T/real(T), T, T, I_N, N, N, Q)

! Q1 matrix
  call KRON(J_T/real(T), T, T, I_N, N, N, QQ)

  do i = 1, nobs
    We(i) = dot_product(ItW(i,:), e)
  end do
  do i = 1, nobs
    WWe(i) = dot_product(ItW(i,:), We)
  end do

RETURN
END SUBROUTINE SETUP_PANELGM



SUBROUTINE KRON(A, p, q, B, r, s, C)
IMPLICIT none
INTEGER, INTENT(IN):: p, q, r, s
REAL, DIMENSION(p,p), INTENT(IN) :: A
REAL, DIMENSION(r,r), INTENT(IN) :: B
REAL, DIMENSION(p*r, p*r), INTENT(OUT) :: C
INTEGER :: h, k, i, j
!**************************************************
! Kronecker product of A(p by q) and B(r by s)
!*************************************************
! same dimension case
!  do k = 1, p
!    do i = 1, r
!      do j = 1, r
!        C(i+(h-1)*r,j+(k-1)*r) = A(h,k)*B(i,j)
!      end do
!    end do
!  end do
!end do

! different dimension case
do h = 1,p
  do k = 1, q
    do i = 1, r
      do j = 1, s
        C(i+(h-1)*r,j+(k-1)*s) = A(h,k)*B(i,j)
      end do
    end do
  end do
end do

RETURN
END SUBROUTINE KRON



SUBROUTINE MOMENTS
USE DataParameters
IMPLICIT none
INTEGER :: i, j
REAL :: con1, con2, con3, con4


 do i = 1, nobs
    Qe(i) = dot_product(Q(i,:), e)
    QWe(i) = dot_product(Q(i, :), We)
    QWWe(i) = dot_product(Q(i,:), WWe)
    QQe(i) = dot_product(QQ(i,:), e)
    QQWe(i) = dot_product(QQ(i,:), We)
    QQWWe(i) = dot_product(QQ(i,:), WWe)
 end do

    eQe = dot_product(e, Qe)
    eQWe = dot_product(e, QWe)
    eQWWe = dot_product(e, QWWe)
    WeQWe = dot_product(We, QWe)
    WWeQWe = dot_product(WWe, QWe)
    WeQWWe = dot_product(We, QWWe)
    WWeQWWe = dot_product(WWe, QWWe)

    eQQe = dot_product(e, QQe)
    eQQWe = dot_product(e, QQWe)
    eQQWWe = dot_product(e, QQWWe)
    WeQQWe = dot_product(We, QQWe)
    WWeQQWe = dot_product(WWe, QQWe)
    WeQQWWe = dot_product(We, QQWWe)
    WWeQQWWe = dot_product(WWe, QQWWe)

  trWTW = 0.0
  do i = 1, N
    do j = 1, N
      trWTW = trWTW + W(i,j)**2
   end do
 end do

! Elements of samples moments
 Gn = 0.0
 con1 = 2.0/ (real(N)*(real(T)-1.0))
 con2 = 1.0/(real(N)*(real(T)-1.0))
 con3 = 2.0/real(N)
 con4 = 1.0/real(N)
 Gn(1,1) = con1*eQWe
 Gn(1,2) = -con2*WeQWe
 Gn(1,3) = 1.0
 Gn(2,1) = con1*WWeQWe
 Gn(2,2) = -con2*WWeQWWe
 Gn(2,3) = con4*trWTW
 Gn(3,1) = con2*(eQWWe + WeQWe)
 Gn(3,2) = -con2*WeQWWe
 Gn(4,1) = con3*eQQWe
 Gn(4,2) = -con4*WeQQWe
 Gn(4,4) = 1.0
 Gn(5,1) = con3*WWeQQWe
 Gn(5,2) = -con4*WWeQQWWe
 Gn(5,3) = con4*trWTW
 Gn(6,1) = con4*(eQQWWe + WeQQWe)
 Gn(6,2) = -con4*WeQQWWe

 small_gn = 0.0
 small_gn(1) = con2*eQe
 small_gn(2) = con2*WeQWe
 small_gn(3) = con2*eQWe
 small_gn(4) = con4*eQQe
 small_gn(5) = con4*WeQQWe
 small_gn(6) = con4*eQQWe

RETURN
END SUBROUTINE MOMENTS



SUBROUTINE NEWTON
USE DataParameters
IMPLICIT none
INTEGER ::       i, iter, j
REAL ::                                    SSE, newSSE

   delta0 = delta_temp

   call SSEfunc(SSE)

 iter = 12

 DO i = 1, iter
   print *, 'iteration', i
   call direction

   delta0 = delta_temp + ddelta
   call SSEfunc(newSSE)
   if (newSSE.lt.SSE) then
     delta_temp = delta_temp + ddelta
     print *, 'newSSE = ', newSSE
     SSE = newSSE
   else
     print *, 'adjusting step length'
     call step(SSE, newSSE)
     delta_temp = delta_temp + ddelta
     SSE = newSSE
   endif

 END DO

RETURN
END SUBROUTINE NEWTON



SUBROUTINE SSEfunc(SSE)
USE DataParameters
IMPLICIT none
INTEGER ::  i, iter
REAL, INTENT(OUT) ::     SSE

! Calculate SSE (Sum of Squared Error)

IF (estimate == 'initialGM') then
   do i = 1,  n_moments0
     Error(i) =Gn0(i,1)*delta0(1) + Gn0(i,2)*delta0(1)**2.0 +&
                          Gn0(i,3)*delta0(2)-small_gn0(i)
     SE(i) = Error(i)**2.0
   end do
 
    SSE = sum(SE)

 ELSE IF ((estimate == 'WeightedGM') .or. (estimate == 'PWeightedGM')) then
    do i = 1, n_moments
     Error(i) =Gn0(i,1)*delta0(1) + Gn0(i,2)*delta0(1)**2.0 +&
              Gn0(i,3)*delta0(2) + Gn0(i,4)*delta0(3)-small_gn0(i)
     SE(i) = Error(i)**2.0
    end do

    SSE = sum(SE)
 END IF 

RETURN
END SUBROUTINE SSEfunc



SUBROUTINE direction
USE DataParameters
IMPLICIT none
INTEGER ::  i

IF (estimate == 'initialGM') then

 ! gradient g
   g = 0.0
   do i = 1, n_moments0 
     g(1) = g(1) + 2.0*Error(i)*(Gn0(i,1)+2.0*delta0(1)*Gn0(i,2))
     g(2) = g(2) + 2.0*Error(i)*Gn0(i,3)
   end do

 ! hessian H
   H = 0.0
   do i = 1, n_moments0 
     H(1,1) = H(1,1) + 2.0*(Gn0(i,1)+2.0*delta0(1)*Gn0(i,2))**2.0 &
                 + 4.0*Error(i)*Gn0(i,2)
     H(1,2) = H(1,2) + 2.0*(Gn0(i,1)+2.0*delta0(1)*Gn0(i,2))*Gn0(i,3)
     H(2,2) = H(2,2) + 2.0*Gn0(i,3)**2.0
   end do
     H(2,1) = H(1,2)

  Hinverse = identity_H
  call la_sysv(H, Hinverse)

! Search direction
  ddelta =-matmul(Hinverse, g)

ELSE IF ((estimate == 'WeightedGM') .or. (estimate == 'PWeightedGM')) then

  ! gradient g
   g = 0.0
   do i = 1, n_moments
     g(1) = g(1) + 2.0*Error(i)*(Gn0(i,1)+2.0*delta0(1)*Gn0(i,2))
     g(2) = g(2) + 2.0*Error(i)*Gn0(i,3)
     g(3) = g(3) + 2.0*Error(i)*Gn0(i,4)
   end do

 ! hessian H
   H = 0.0
   do i = 1, n_moments
     H(1,1) = H(1,1) + 2.0*(Gn0(i,1)+2.0*delta0(1)*Gn0(i,2))**2.0 &
                 + 4.0*Error(i)*Gn0(i,2)
     H(1,2) = H(1,2) + 2.0*(Gn0(i,1)+2.0*delta0(1)*Gn0(i,2))*Gn0(i,3)
     H(1,3) = H(1,3) + 2.0*(Gn0(i,1)+2.0*delta0(1)*Gn0(i,2))*Gn0(i,4)
     H(2,2) = H(2,2) + 2.0*Gn0(i,3)**2.0
     H(2,3) = H(2,3) + 2.0*Gn0(i,4)*Gn0(i,3)
     H(3,3) = H(3,3) + 2.0*Gn0(i,4)**2.0
   end do
     H(2,1) = H(1,2)
     H(3,1) = H(1,3)
     H(3,2) = H(2,3)

  Hinverse = identity_H
  call la_sysv(H, Hinverse)

  ddelta =-matmul(Hinverse, g)

END IF

RETURN
END SUBROUTINE direction



SUBROUTINE step(SSE, newSSE)
USE  DataParameters 
IMPLICIT none
REAL, INTENT(IN)    :: SSE
REAL, INTENT(INOUT) :: newSSE

DO WHILE(newSSE < SSE-1d-8)
   print *, 'SSE and last newSSE', SSE, newSSE

   ddelta = 0.8d0*ddelta
   delta0 = delta_temp + ddelta
   call SSEfunc(newSSE)
END DO

RETURN
END SUBROUTINE step



SUBROUTINE  WEIGHTED_MOMENTS
USE DataParameters
IMPLICIT none
INTEGER ::    i, j

! calculate Tw matrix
IF (estimate == 'WeightedGM') then
 ! calculate traces
 ! Tw(2,2)
   do i = 1, N
    do j = 1, N
      WTW(i,j) = dot_product(W(:,i),W(:,j))
    enddo
   enddo

   do i = 1, N
    do j = 1, N
     WTWWTW(i,j) = dot_product(WTW(i,:), WTW(:,j))
    end do
   end do

   trWTWWTW = 0.0
   do i = 1, N 
    trWTWWTW = trWTWWTW + WTWWTW(i,i)
   end do

 ! Tw(2,3) and Tw(3,2)
   do i = 1, N
    do j = 1, N
     WTplusW(i,j) = W(j,i) + W(i,j)
    end do
  end do

  do i = 1, N
   do j = 1, N
    WTWWTplusW(i,j) = dot_product(WTW(i,:), WTplusW(:,j))
   end do
  end do

  trWTWWTplusW = 0.0
  do i = 1, N 
   trWTWWTplusW = trWTWWTplusW + WTWWTplusW(i,i)
  end do

 ! Tw(3,3)
  trWW = 0.0
  do i = 1, N
   do j = 1, N
     trWW = trWW + W(i,j)*W(j,i)
   end do
  end do
 trWWplusWTW = trWW + trWTW

! Tw matrix
   Tw(1,1) = 2
   Tw(1,2) = (2.0/N)*trWTW
   Tw(2,1) = Tw(1,2)
   Tw(2,2) = (2.0/N)*trWTWWTW
   Tw(2,3) = trWTWWTplusW/N
   Tw(3,1) = 0.0
   Tw(3,2) = Tw(2,3)
   Tw(3,3) = trWWplusWTW/N

ELSE IF (estimate == 'PWeightedGM') then

   Tw = 0.0
   do i = 1, n_param
     Tw(i,i) = 1.0
   end do

END IF

 call KRON(Xitemp, 2, 2, Tw, 3, 3, Xi)
 XiInverse = identity_Xi
 call la_posv(Xi, XiInverse)

 ! cholesky decomposition of XiInverse
 call la_potrf(XiInverse, UPLO = 'L')
 do i = 1, n_moments-1
   do j = i+1, n_moments
    XiInverse(i, j) = 0.0
   end do
 end do

 ! Transformation of Gn and small_gn
   Gntemp = Gn
   Gn = 0.0
   Gn = matmul(XiInverse, Gntemp)

   small_gn_temp = small_gn
   small_gn = 0.0
   small_gn = matmul(XiInverse, small_gn_temp)

RETURN
END SUBROUTINE WEIGHTED_MOMENTS



SUBROUTINE GLS
USE DataParameters
IMPLICIT none
INTEGER::          i, j
REAL(wp), DIMENSION(KVars,KVars):: XX
!*****************************************************************
! Solve for GLS estimator of beta, using final rho value
!*****************************************************************
! First transformation of Y and X
 call KRON(I_T, T, T, I_N - rho*W , N, N, IIrhoW)
 do i = 1, nobs 
  do j = 1, KVars
   XTilde(i,j) = dot_product(IIrhoW(i,:), X(:,j))
  end do
  YTilde(i) = dot_product(IIrhoW(i,:), Y)
 end do

! Second transformation of Y and X
 OmegInv = (1.0/sqrt(sigma_sq))*Q + (1.0/sqrt(sigma_one_sq))*QQ
 do i = 1, nobs
   do j = 1, KVars 
     XXTilde(i,j) = dot_product(OmegInv(i,:), XTilde(:,j))
  end do
  YYTilde(i) = dot_product(OmegInv(i,:), YTilde)
 end do

 do i = 1, KVars
   do j = 1, KVars
     XTX(i,j) = dot_product(XXTilde(:,i),XXTilde(:,j))
   enddo
   XTY(i) = dot_product(XXTilde(:,i),YYTilde)
 enddo

  XX = XTX

  call la_posv(XTX,XTY)
  beta = XTY
  print *, "Final beta is ", beta

  etilde = YYTilde - matmul(XXTilde,beta)
  sigma_sq = dot_product(etilde,etilde)/real(nobs)
  print *, "Final sigma squared is ", sigma_sq

  ! Invert XTX
  call la_posv(XX,identity_beta)

  ! Multiply the inverse by sigma-squared to obtain
  !  the variance matrix for beta
  var_beta = sigma_sq*identity_beta

  do i = 1, KVars
   print *, beta(i), sqrt(var_beta(i,i)), beta(i)/sqrt(var_beta(i,i))
  enddo


RETURN
END SUBROUTINE GLS


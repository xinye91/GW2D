
subroutine petsc2DSolver(m,n,M5,RHS,tol,hNew)
   use proc_mod
   implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscviewer.h>

   ! implicit none

   ! input argument
   integer :: m, n ! The domain we solved is m*n, and the matrix would be (m*n)x(m*n) 
   REAL*8,DIMENSION(m*n,5) :: M5 ! The five diagonal matrices
   REAL*8,DIMENSION(m*n) :: RHS ! The RHS vector
   REAL*8,DIMENSION(:), pointer :: hNew ! result vector on the processor part
   real*8 tol

   PetscErrorCode ierr
   Mat A
   Vec x,b
   KSP ksp
   PetscInt istart, iend, ione, nrow
   PetscInt ii, jj, i, j, its
   PetscMPIInt rank, nproc   

!   Integer :: md, dm, begi
   
   ione = 1

   ! petsc env setting
   call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
   call MPI_Comm_size(PETSC_COMM_WORLD,nproc,ierr)

   ! row decomposition
   ! now it is handled by proc_mod
!   md = mod(n*m,nproc)
!   dm = n*m/nproc
!   if (rank < md) then
!      begi = rank * (dm+1) + 1
!      nrow = dm + 1
!   else
!      begi = md * (dm+1) + dm*(rank-md) + 1
!      nrow = dm
!   endif
   ! call proc_mod
   nrow = pci%nrow(rank+1)   

   write(*,*)rank,  nrow
 

   ! construct petsc Mat
   call MatCreate(PETSC_COMM_WORLD,A,ierr)
   call MatSetSizes(A,nrow,nrow,m*n,m*n,ierr)
   call MatSetFromOptions(A,ierr)
   call MatSetUp(A,ierr)

   call MatGetOwnershipRange(A,Istart,Iend,ierr)
   write(*,*) rank, ':', Istart, Iend
   allocate(hNew(nrow))   

   do ii=istart,iend-1
      i = II/m
      j = II - i*m
      if (i.gt.0) then
         JJ = II - m
         call MatSetValues(A,ione,II,ione,JJ,M5(ii+1,1),INSERT_VALUES,ierr)
      endif
      if (i.lt.n-1) then
         JJ = II + m
         call MatSetValues(A,ione,II,ione,JJ,M5(ii+1,5),INSERT_VALUES,ierr)
      endif
      if (j.gt.0) then
         JJ = II - 1
         call MatSetValues(A,ione,II,ione,JJ,M5(ii+1,2),INSERT_VALUES,ierr)
      endif
      if (j.lt.m-1) then
         JJ = II + 1
         call MatSetValues(A,ione,II,ione,JJ,M5(ii+1,4),INSERT_VALUES,ierr)
      endif
      call MatSetValues(A,ione,ii,ione,ii,M5(ii+1,3),INSERT_VALUES,ierr)
   end do
write(*,*) "?????"
   call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
write(*,*) '!!!!!'
   call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
write(*,*) '.......'

   !call Matview(A,PETSC_VIEWER_STDOUT_WORLD,ierr) 

   ! construct petsc vector
   call VecCreateMPI(PETSC_COMM_WORLD,nrow,m*n,b,ierr)
   call VecSetFromOptions(b,ierr)


   call VecGetOwnershipRange(b,istart,iend, ierr)
write(*,*) Istart, Iend

   do ii=istart,iend-1
      call VecSetValues(b,ione,ii,RHS(ii+1),INSERT_VALUES,ierr) 
   enddo 
   call VecAssemblyBegin(b,ierr)
   call VecAssemblyEnd(b,ierr)  
   call VecDuplicate(b,x,ierr) 

  ! call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr) 

istart = 100
   call VecGetOwnershipRange(x,istart,iend, ierr)
write(*,*) Istart, Iend

   ! constrcut KSP
   call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
   call KSPSetOperators(ksp,A,A,ierr)
!   call KSPGetPC(ksp,pc,ierr)
!   ptype = PCJACOBI
!   call PCSetType(pc,ptype,ierr)
  
   call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_REAL,     &
     &     PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)  
   call KSPSetFromOptions(ksp,ierr)
   
   ! solve the linear system
   call KSPSolve(ksp,b,x,ierr)
  
   call KSPGetIterationNumber(ksp,its,ierr)
   if (rank .eq. 0) then
      write(*,*) 'iterations', its
   endif
  
   ! call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr) 
   ! copy to real data
   do ii = istart, iend-1 
     call VecGetValues(x,ione,ii,hNew(ii+1-istart),ierr)
   enddo
   write(*,*) rank, 'VecGet:', hNew

   call KSPDestroy(ksp,ierr)
   call VecDestroy(x,ierr)
   call VecDestroy(b,ierr)
   call MatDestroy(A,ierr)


end subroutine


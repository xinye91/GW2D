
module proc_mod

   use gd_mod
   include 'mpif.h'
   private

   public :: proc_info, pci
   public :: createPCInfo_2D, linkPciFromGd

   type proc_info
     integer, pointer :: nx(:), nrow(:)
     integer, pointer :: Xrange(:,:), Yrange(:,:)
   end type proc_info

   type(proc_info), save, target :: pci

  contains

   subroutine linkPciFromGd(gd_rt) 
      implicit none
      type(gd), target :: gd_rt
      type(gd), pointer :: gdT
      
      gdT => gd_rt .G. 'proc'

      call getPtr(gdT,'nx',pci%nx)
      call getPtr(gdT,'nrow',pci%nrow)
      call getPtr(gdT,'Yrange',pci%Yrange)
      call getPtr(gdT,'Xrange',pci%Xrange)
      

   end subroutine

   subroutine createPCInfo_2D(domainSize)

      implicit none
      integer :: domainSize(2)
      integer :: nproc, rank, md, dm
      integer :: i, ierr
      integer, pointer :: nx(:), nrow(:), Yrange(:,:), Xrange(:,:)

      ! 2D domain decomposition
      ! fortran stores 2nd dim consecutively.
      ! we regard 1st dim as y, 2nd dim as x.
      call MPI_Init(ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
      allocate(nx(0:nproc-1))
      allocate(nrow(0:nproc-1))
      allocate(Yrange(2,0:nproc-1))
      allocate(Xrange(2,0:nproc-1))
      md = mod(domainSize(2)-2,nproc)
      dm = (domainSize(2)-2)/nproc
      do i=0,nproc-1
         ! nx: the number of columns the proc gets (not including ghost cells)
         if (i < md) then
            nx(i) = dm + 1
         else
            nx(i) = dm
         endif
         ! Yrange, Xrange: including ghost cell, the range of Y(1st), X(2nd) in ths proc's domain
         Yrange(1,i)=1; Yrange(2,i)=domainSize(1)  ! including ghost cells: 1 - msize(1). same for all procs.
         if (i==0) then
            Xrange(1,i)=1; Xrange(2,i)=nx(i)+2 
            ! including ghost cells: 1 - nx(0)+2, where:
            ! 1 is the domain's ghost cell in this direction,
            ! nx(0)+1 is 'calculating domain''s ending, so nx(0)+2 is the ghost cell for proc #0
         else
            Xrange(1,i)=Xrange(2,i-1)-1; Xrange(2,i)=Xrange(1,i)+nx(i)+1
            ! The left ghost cell is the last proc's last 'calculating end'
            ! and the decompoced domain length in this direction is nx(i)+2
         endif
         nrow(i)=(domainSize(1)-2)*nx(i)
      enddo
      ! put into gd
      i = gd_base%np + 1
      gd_base%g(i)%name = 'proc'
      gd_base%f(i) = 'proc'
      call allocate_gd(gd_base%g(i))
      gd_base%np = i
      !write(*,*) i
      call addptr((/nproc/),'nx',nx,gd_base%g(i))   
      call addptr((/nproc/),'nrow',nrow,gd_base%g(i))   
      call addptr((/2,nproc/),'Yrange',Yrange,gd_base%g(i))   
      call addptr((/2,nproc/),'Xrange',Xrange,gd_base%g(i))

   end subroutine

end module proc_mod



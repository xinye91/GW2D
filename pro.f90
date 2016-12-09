program Ncgd

   use netcdf
   use gd_mod
   use displayMod
   use ncTableHier
   use vdata
   use flow
   use proc_mod
   implicit none
   include 'mpif.h'



   type(gd) :: gdT, gdT2, gdT3
   type(gd), pointer :: gr_ptr
   real*8, pointer :: msize(:,:)
   type(VDZ_type), pointer :: VDZ
   real*8, pointer :: MAPP(:,:)
   integer :: s(2), nGh(2)
   integer :: NBC

   call buildNCTable('example2.nc')
   write(*,*) '==='   
   gdT = gd_base .G. 'example2.nc'

   gdT2 = gdT .G. 'g'
   !gdT3 = gdT2 .G. 'DM'  
   call getPtr(gdT2,'GW',gr_ptr,2)
   call getPtr(gr_ptr,'K',msize)
!   write(*,*) msize
   call getPtr(gdT2,'GW',gr_ptr,1)
   call getPtr(gr_ptr,'K',msize)
!   write(*,*) msize

   ! link
   call linkVdataFromGd(gdT)   
   !write(*,*) w%DM%mask
   s(1)=size(w%DM%mask,1) 
   s(2)=size(w%DM%mask,2)

   write(*,*) s
   call createPCInfo_2D(s)   
   write(*,*)  'proc'

   call linkPciFromGd(gd_base)
   !write(*,*) pci%Xrange(:,2)

   nGh = (/1,1/)
   allocate(MAPP(s(1),s(2)))
   MAPP = 0.0D0
   where(W%DM%mask)
     MAPP = 1.0D0
   end where
   !call Aquifer(g%GW(1),MAPP,g%GW(1)%Cond,size(g%GW(1)%Cond),1)
   NBC=size(g%GW(1)%Cond)
   call Aquifer2D(g%GW(1)%ny,g%GW(1)%nx,nGh,g%GW(1)%h,g%GW(1)%EB,g%GW(1)%D,g%GW(1)%K,g%GW(1)%W,      &
     &          g%GW(1)%ST, Mapp, g%GW(1)%Cond,NBC,g%GW(1)%dd,g%GW(1)%dt,g%GW(1)%hNew,g%GW(1)%DR, 1) 
   
!   write(*,*) g%GW(1)%h

end program

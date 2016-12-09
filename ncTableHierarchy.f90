MODULE ncTableHier
    ! netcdf file interface
    ! modified form matTable
    ! This code replaces matTable.f90 when we link to netcdf (in serial) instead of .mat file
    ! *****************************************************
    ! Xinye Ji
    ! Pennsylvania State University
    ! Civil and Environmental Engineering
    ! jixinye91@gmail.com
    ! *****************************************************
    ! build an intermediate table object which stores data read from a netcdf file. 
    ! It can be used later to conveniently link the data to computational objects
    ! syntax:
    ! ####################################################
    ! SOME WHERE IN THE PROGRAM YOU MUST FIRST DO:
    ! ++ use ncTableHier
    ! ++ call buildNCTable(nc_file)
    ! ####################################################
    ! after this, you can access all variable stored in the nc_file via gd_base, which is a type(gd) variable in module ncTable
    ! the general syntax is
    ! ++ call getPtr(gd_base,fieldname,ptr,num)  ! num is optional, default to 1
    ! this will link the data to the pointer, ptr, which may be a real*8,integer, logical array supported up to 6 dimension
    ! or a type(gd) array pointer.
    ! depending on what pointer is passed in, the calling interface will pick the corresponding type/dimension subroutine
    ! If the queried field is not there or if the field type/dimension does not match that of the ptr, an error is thrown.
    ! num is the index, if omitted, default to 1. type(gd) allows to access sub-structures.
    ! Overloaded operators as a fast access:
    ! Besides the calling interface getPtr, we have overloaded operator .G., .FF. and .FFF. serving as shorthand methods
    ! for extracting sub-structure, [real*8,dimension(:,:)] and [real*8,dimension(:,:,:)] data, respectively. 
    ! With these operators, num is always 1.
      ! *****************************************************
      !Copyright (c) 2015, Chaopeng Shen, Xinye Ji
      !All rights reserved.
      !
      !THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
      !ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
      !WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
      !DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
      !ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
      !(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
      !LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
      !ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
      !(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
      !SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
      !    
      ! Please cite the following two papers if this tool is used in your academic research work
      ! 1. Ji, XY and CP. Shen, The introspective may achieve more: enhancing existing Geoscientific 
      !      models with native-language structural reflection, Environmental Modeling & Software.
      ! 2. Shen, CP., J. Niu and K. Fang, Quantifying the Effects of Data Integration Algorithms
      !      on the Outcomes of a Subsurface - Land Surface Processes Model, 
      !      Environmental Modeling & Software., doi: 10.1016/j.envsoft.2014.05.006 (2014)
      !
      !Commercial use of the software need to obtain explicit consent from the copyright owner.
      !Redistribution and use in source and binary forms for non-commercial purposes, with or without
      !modification, are permitted provided that the following conditions are met:
      !
      !1. Redistributions of source code must retain the above copyright notice, this
      !   list of conditions and the following disclaimer.
      !2. Redistributions in binary form must reproduce the above copyright notice,
      !   this list of conditions and the following disclaimer in the documentation
      !   and/or other materials provided with the distribution.
      !
      !*******************************************************      
    use gd_mod
    use netcdf
    
    integer nG, nsG
    parameter (nG =10)
    parameter (nsG = 60)
    integer, parameter::maxStrLen=100
    PUBLIC :: buildNCTable, connectNCField
    contains
    
    subroutine buildNCTable(nc_file)
      ! read a workspace saved in nc_file into (gd_base .G. nc_file)
      implicit none
      integer, parameter :: cLen = 100
      integer*4 num, nbase
      character*(cLen) names(6000), name
      integer   i, status, np
      integer ntemp,j
      character*(*) nc_file
      character*(cLen), dimension(10) :: skipVar
      integer nSkipVar, skip
      integer :: ncid, varid, dimid, numVar


      nSkipVar = 1
      skipVar(1) = 'Prob'

      status = nf90_open(nc_file, NF90_NOWRITE, ncid)
      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          stop "Stopped in nf90_open"
      end if

      status = nf90_inquire(ncid, nvariables = numVar)
      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          stop "Stopped in nf90_inquire"
      end if

      ! Get variable name
      do i=1,numVar
          status = nf90_inquire_variable(ncid,i,name = names(i))
          if (status /= nf90_noerr) then 
              print *, trim(nf90_strerror(status))
              stop "Stopped in nf90_inquire_variable"
          end if
      enddo
      
      nbase = 10
      if (.not. associated(gd_base%p)) call allocate_gd(gd_base,nbase)
      gd_base%name(1) = 'base'
      np = gd_base%np + 1   
      gd_base%f(np) = nc_file
      gd_base%np    = np
      
      if (.not. associated(gd_base%g)) then
         allocate(gd_base%g(nbase))
      endif
      call allocate_gd(gd_base%g(np),nsG)
      gd_base%g(np)%name = nc_file
      
      do i=1,numVar
         skip = 0
         do j=1, nSkipVar
           if (names(i) .eq. skipVar(j)) skip=1
         enddo
         if (skip .ne. 1) call connectNCField(gd_base%g(np),i,names(i),ncid,nc_file)
      enddo
      
      status = nf90_close(ncid) 
      if (status /= nf90_noerr) then 
          print *, trim(nf90_strerror(status))
          stop "Stopped in nf90_close"
      end if
      
    end subroutine
    
    
    
    
    recursive subroutine connectNCField(main,dat,field,ncid,gdname)
      ! add "field" to the gd table "main" by connecting to "dat"
      ! so that main%f(np+1)=field and main%p(np+1)%ptr => dat (unwrapped from mx)
      ! XY: now dat = varid
      implicit none
      type(gd),target:: main
      integer m,n,j
      character(*),intent(in)::field
      character(*),intent(in),optional::gdname
      integer :: dat, ncid
      
      integer*4 fieldnum,nfields,nallo
      integer*4 i
      integer :: ndims, status, varid, xtype
      integer :: dimids(NF90_MAX_VAR_DIMS)
      logical :: islogic, isreal, isint, ischar
      
      real*8, pointer :: data_R1(:)
      real*8, pointer :: data_R2(:,:)
      real*8, pointer :: data_R3(:,:,:)
      real*8, pointer :: data_R4(:,:,:,:)
      real*8, pointer :: data_R5(:,:,:,:,:)
      real*8, pointer :: data_R6(:,:,:,:,:,:)
      integer, pointer :: data_I1(:)
      integer, pointer :: data_I2(:,:)
      integer, pointer :: data_I3(:,:,:)
      integer, pointer :: data_I4(:,:,:,:)     
      integer, pointer :: data_I5(:,:,:,:,:)
      integer, pointer :: data_I6(:,:,:,:,:,:)    
      
      integer di(6),ns(6)
      integer numel
      integer :: dims(6) = (/0, 0, 0, 0, 0, 0/)
      integer indD, indB1, indB2, Seq
      character(maxStrLen) :: subname, mdName, mdName2, mdSeq
      character*(maxStrLen), pointer :: namestr
      
      
      if (dat .eq. 0) return ! empty!!

      varid = dat
      isreal  = .false.
      islogic = .false.
      isint   = .false.
      ischar  = .false.
      
      status = nf90_inquire_variable(ncid, varid, xtype=xtype, ndims=ndims, dimids=dimids)
      
      do i = 1, ndims
          status = nf90_inquire_dimension(ncid, dimids(i), len=dims(i))
      end do

      numel = 1; ns = 0
      DO j=1,ndims
        numel = numel*dims(j)
        if (dims(j) > 1) ns(j)=1
      ENDDO
      
      di = dims
      
      if (numel .eq. 0) return ! empty!
      
      ! Hierarchy judgement
      indD = index(field,'.')
      if (indD > 0) then! there is hierarchy
          mdname = field(1:indD-1)
          subname = field(indD+1:)
        
          allocate(namestr)
          if (present(gdname)) then
              namestr = trim(gdname) // trim('.')// trim(mdname)
          else
              namestr = trim(main%name(1)) // trim('.')// trim(mdname)
          endif 
          
          ! structure array judgement
          indB1 = index(mdname,'(')
          Seq = 1
          if (indB1 > 0) then
              indB2 = index(mdname,')')
              mdname2 = mdname(1:indB1-1)
              mdSeq = mdname(indB1+1:indB2-1)
              Seq = str2num(trim(mdSeq))
              mdname = mdname2
          end if
          
          ! search if it exists
          i = searchStr(main%f,mdname,main%np,Seq)
          if (i<=0) then ! create it
              if (.not. associated(main%g)) then
                  allocate(main%g(nsG))
              endif
              n = main%np + 1
              call allocate_gd(main%g(n))
              main%g(n)%name(1) = namestr
              main%f(n) = mdname
              main%np = n
          else ! point to it
              n = i
          endif
          call connectNCField(main%g(n), dat, subname, ncid)
              
      else
          ! datatype and get data
          if ((xtype .eq. NF90_DOUBLE) .or. (xtype .eq. NF90_FLOAT)) then
              isreal = .true.
              n = main%np + 1
              select case (ndims)
              case (1)
                  allocate(data_R1(dims(1)))
                  status = nf90_get_var(ncid, varid, data_R1)              
                  if (status /= nf90_noerr) then                   
                      print *, trim(nf90_strerror(status))                  
                      stop "Stopped in nf90_get_var"
                  end if
                  call assignR1(numel,data_R1,main,n)
              case (2)
                  allocate(data_R2(dims(1),dims(2)))
                  status = nf90_get_var(ncid, varid, data_R2)
                  if (status /= nf90_noerr) then 
                      print *, trim(nf90_strerror(status))
                      stop "Stopped in nf90_get_var"
                  end if
                  call assignR2(di,data_R2,main,n)
              case (3)
                  allocate(data_R3(dims(1),dims(2),dims(3)))
                  status = nf90_get_var(ncid, varid, data_R3)
                  if (status /= nf90_noerr) then 
                      print *, trim(nf90_strerror(status))
                      stop "Stopped in nf90_get_var"
                  end if
                  call assignR3(di,data_R3,main,n)
              case (4)
                  allocate(data_R4(dims(1),dims(2),dims(3),dims(4)))
                  status = nf90_get_var(ncid, varid, data_R4)
                  if (status /= nf90_noerr) then 
                      print *, trim(nf90_strerror(status))
                      stop "Stopped in nf90_get_var"
                  end if
                  call assignR4(di,data_R4,main,n)
              case (5)
                  allocate(data_R5(dims(1),dims(2),dims(3),dims(4),dims(5)))
                  status = nf90_get_var(ncid, varid, data_R5)
                  if (status /= nf90_noerr) then 
                      print *, trim(nf90_strerror(status))
                      stop "Stopped in nf90_get_var"
                  end if
                  call assignR5(di,data_R5,main,n)
              case (6)
                  allocate(data_R6(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
                  status = nf90_get_var(ncid, varid, data_R6)
                  if (status /= nf90_noerr) then 
                      print *, trim(nf90_strerror(status))
                      stop "Stopped in nf90_get_var"
                  end if
                  call assignR6(di,data_R6,main,n)  
              end select
              main%f(n) = field
              main%np = n
          elseif ((xtype .eq. NF90_SHORT) .or. (xtype .eq. NF90_INT) .or.  (xtype .eq. NF90_BYTE)) then
              isint = .true.
              n = main%np + 1
              select case (ndims)
              case (1)
                  allocate(data_I1(dims(1)))
                  status = nf90_get_var(ncid, varid, data_I1)
                  if (status /= nf90_noerr) then 
                      print *, trim(nf90_strerror(status))
                      stop "Stopped in nf90_get_var"
                  end if
                  call assignI1(numel,data_I1,main,n)
              case (2)
                  allocate(data_I2(dims(1),dims(2)))
                  status = nf90_get_var(ncid, varid, data_I2)
                  if (status /= nf90_noerr) then 
                      print *, trim(nf90_strerror(status))
                      stop "Stopped in nf90_get_var"
                  end if
                  call assignI2(di,data_I2,main,n)
              case (3)
                  allocate(data_I3(dims(1),dims(2),dims(3)))
                  status = nf90_get_var(ncid, varid, data_I3)
                  if (status /= nf90_noerr) then 
                      print *, trim(nf90_strerror(status))
                      stop "Stopped in nf90_get_var"
                  end if
                  call assignI3(di,data_I3,main,n)
              case (4)
                  allocate(data_I4(dims(1),dims(2),dims(3),dims(4)))
                  status = nf90_get_var(ncid, varid, data_I4)
                  if (status /= nf90_noerr) then 
                      print *, trim(nf90_strerror(status))
                      stop "Stopped in nf90_get_var"
                  end if
                  call assignI4(di,data_I4,main,n)
              case (5)
                  allocate(data_I5(dims(1),dims(2),dims(3),dims(4),dims(5)))
                  status = nf90_get_var(ncid, varid, data_I5)
                  if (status /= nf90_noerr) then 
                      print *, trim(nf90_strerror(status))
                      stop "Stopped in nf90_get_var"
                  end if
                  call assignI5(di,data_I5,main,n)
              case (6)
                  allocate(data_I6(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
                  status = nf90_get_var(ncid, varid, data_I6)
                  if (status /= nf90_noerr) then 
                      print *, trim(nf90_strerror(status))
                      stop "Stopped in nf90_get_var"
                  end if
                  call assignI6(di,data_I6,main,n)  
              end select
              main%f(n) = field
              main%np = n
          elseif (xtype .eq. NF90_CHAR) then
              ischar = .true.
              ! temporary no char data.
          end if
      endif ! hierarchy judgement
      
    end subroutine

    
END MODULE ncTableHier

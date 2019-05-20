!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  subroutine write_output_ASCII_or_binary(one_seismogram, &
                                          NSTEP,it,SIMULATION_TYPE,DT,t0, &
                                          iorientation,sisname,final_LOCAL_PATH)

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function

  use constants

  use specfem_par, only: myrank,USE_BINARY_FOR_SEISMOGRAMS,SAVE_ALL_SEISMOS_IN_ONE_FILE,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
    seismo_offset,seismo_current

  implicit none

  integer,intent(in) :: NSTEP,it,SIMULATION_TYPE

  real(kind=CUSTOM_REAL), dimension(NDIM,NTSTEP_BETWEEN_OUTPUT_SEISMOS),intent(in) :: one_seismogram

  double precision,intent(in) :: t0,DT

  integer,intent(in) :: iorientation

  character(len=MAX_STRING_LEN),intent(in) :: sisname,final_LOCAL_PATH

  ! local parameter
  integer :: isample,ier,it_current
  real(kind=CUSTOM_REAL) :: value,time_t

  real, dimension(1:seismo_current) :: tr

  ! opens seismogram file
  if (USE_BINARY_FOR_SEISMOGRAMS) then

    ! binary format case
    if (SAVE_ALL_SEISMOS_IN_ONE_FILE) then
      ! trace name
      write(IOUT) sisname
    else
      if (seismo_offset == 0) then
        open(unit=IOUT, file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)), &
             form='unformatted', action='write', access='stream', status='replace', iostat=ier)
      else
        open(unit=IOUT, file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)), &
             form='unformatted', action='write', access='stream', status='old', position='append', iostat=ier)
      endif
      if (ier /= 0) call exit_mpi(myrank,'Error opening file: '// &
                                  final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)))
    endif

  else

    ! ASCII format case
    if (SAVE_ALL_SEISMOS_IN_ONE_FILE) then
      ! trace name
      write(IOUT,*) sisname(1:len_trim(sisname))
    else
      if (seismo_offset == 0) then
        open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)), &
             action='write', status='unknown',iostat=ier)
      else
        open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)), &
             action='write', status='old',position='append',iostat=ier)
      endif
      if (ier /= 0) call exit_mpi(myrank,'Error opening file: '// &
                                  final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)))
    endif

  endif

  ! make sure we never write more than the maximum number of time steps
  ! subtract half duration of the source to make sure travel time is correct
  do isample = 1,seismo_current

    ! seismogram value
    value = one_seismogram(iorientation,isample)

    ! current time increment
    it_current = seismo_offset + isample
    if (it_current > it) stop 'Invalid current time step in writing output ASCII/binary'

    ! time
    if (SIMULATION_TYPE == 1) then
      ! forward simulation
      ! distinguish between single and double precision for reals
      time_t = real( dble(it_current-1)*DT - t0 ,kind=CUSTOM_REAL)
    else if (SIMULATION_TYPE == 3) then
      ! adjoint simulation: backward/reconstructed wavefields
      ! distinguish between single and double precision for reals
      ! note: compare time_t with time used for source term
      time_t = real( dble(NSTEP-it_current)*DT - t0 ,kind=CUSTOM_REAL)
    endif

    if (USE_BINARY_FOR_SEISMOGRAMS) then
      ! binary format case
      if (SAVE_ALL_SEISMOS_IN_ONE_FILE) then
        write(IOUT) time_t,value
      else
        tr(isample) = value
      endif
    else
      ! ASCII format case
      write(IOUT,*) time_t,value
    endif

  enddo

  if (.not. SAVE_ALL_SEISMOS_IN_ONE_FILE) then
    ! binary format case
    ! writes out whole trace into binary file
    if (USE_BINARY_FOR_SEISMOGRAMS) write(IOUT) tr(:)
    ! closes files
    close(IOUT)
  endif

  end subroutine write_output_ASCII_or_binary

  subroutine write_output_ASCII_or_binary_der(one_seismogram, &
                                          NSTEP,it,SIMULATION_TYPE,DT,t0, &
                                          iorientation,sisname,final_LOCAL_PATH)

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function

  use constants

  use specfem_par,only: USE_BINARY_FOR_SEISMOGRAMS,SAVE_ALL_SEISMOS_IN_ONE_FILE

  implicit none

  integer :: NSTEP,it,SIMULATION_TYPE
 
  !Carene had to rewrite an entire subroutine because of the next argument
  !passing alone
  real(kind=CUSTOM_REAL), dimension(NDIM+1,NSTEP) :: one_seismogram

  double precision t0,DT

  integer iorientation

  character(len=MAX_STRING_LEN) :: sisname,final_LOCAL_PATH

  ! local parameter
  integer isample,nt_s,ier
  real, dimension(:), allocatable :: tr
  real(kind=CUSTOM_REAL) :: time_t

  ! opens seismogram file
  if (USE_BINARY_FOR_SEISMOGRAMS) then

    ! allocate trace
    nt_s = min(it,NSTEP)
    allocate(tr(nt_s),stat=ier)
    if (ier /= 0) stop 'error allocating array tr()'

    ! binary format case
    if(.not. SAVE_ALL_SEISMOS_IN_ONE_FILE) then
      open(unit=IOUT, file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//&
         sisname(1:len_trim(sisname)), form='unformatted', access='direct', recl=4*(nt_s))
    else
      write(IOUT) sisname
    endif
    tr(:)=0

  else

    ! ASCII format case
    if(.not. SAVE_ALL_SEISMOS_IN_ONE_FILE) then
      open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//&
         sisname(1:len_trim(sisname)),status='unknown')
    else
      write(IOUT,*) sisname(1:len_trim(sisname))
    endif

  endif

  ! make sure we never write more than the maximum number of time steps
  ! subtract half duration of the source to make sure travel time is correct
  do isample = 1,min(it,NSTEP)

      ! forward simulation
      if (SIMULATION_TYPE == 1) then
        ! distinguish between single and double precision for reals
        if (CUSTOM_REAL == SIZE_REAL) then
          time_t = sngl( dble(isample-1)*DT - t0 )
        else
          time_t = dble(isample-1)*DT - t0
        endif
      endif

      ! adjoint simulation: backward/reconstructed wavefields
      if (SIMULATION_TYPE == 3) then
        ! distinguish between single and double precision for reals
        ! note: compare time_t with time used for source term
        if (CUSTOM_REAL == SIZE_REAL) then
          time_t = sngl( dble(NSTEP-isample)*DT - t0 )
        else
          time_t = dble(NSTEP-isample)*DT - t0
        endif
      endif

      if (USE_BINARY_FOR_SEISMOGRAMS) then
        ! binary format case
        if(.not. SAVE_ALL_SEISMOS_IN_ONE_FILE) then
          tr(isample) = one_seismogram(iorientation,isample)
        else
          write(IOUT) time_t,one_seismogram(iorientation,isample)
        endif
      else
        ! ASCII format case
        write(IOUT,*) time_t,' ',one_seismogram(iorientation,isample)
      endif

  enddo

  ! binary format case
  if (USE_BINARY_FOR_SEISMOGRAMS) then
    ! writes out whole trace into binary file
    if(.not. SAVE_ALL_SEISMOS_IN_ONE_FILE) write(IOUT,rec=1) tr
    deallocate(tr)
  endif

  if(.not. SAVE_ALL_SEISMOS_IN_ONE_FILE) close(IOUT)

  end subroutine write_output_ASCII_or_binary_der


! This code generates posX files containing center-of-mass-weighted (COM) coordinates for each time step from the NVE trajectory calculated using VASP, where X is the Xth timestep of the MD simulation.
! Works using VASP 5.3.3
! Notes: Currently coded for systems with three atoms - must modify to extend to more!
! Input: 1) 'input' file: contains a row for each atom type listing atom name as it should appear in Gaussian09 input file, number of atoms of this type, and its atomic mass in atomic units.  The final row should contain the total number of time steps in the MD trajectory and the total number of atoms.  All input should be separated by at least a space.
!        2) 'MD_G09' file: OUTCAR file from VASP NVE simulation.
program readVASP
implicit none

character(len=128) :: t,filename,tmp,atom1,atom2,atom3,COM
integer :: nat, istep, k, i,natom1,natom2,natom3
real, allocatable :: coord(:,:,:)
real, allocatable :: coord_sh(:,:,:)
real, allocatable :: COM_x(:), COM_y(:), COM_z(:)
real :: mass_tot, mass1, mass2, mass3

open(6,file='input_cfg_gen',status='old')
open(7,file='OUTCAR',status='old')

! Read input file with atom names (atom1,atom2,atom3), number (natom1,natom2,natom3), mass (mass1,mass2,mass3), and MD total time steps (istep) and total atom number (nat)
read(6,*) atom1,natom1,mass1
read(6,*) atom2,natom2,mass2
read(6,*) atom3,natom3,mass3
read(6,*) istep,nat

close(6)

! Calculate total mass of system
mass_tot = natom1*mass1 + natom2*mass2 + natom3*mass3

! Allocate matrices
allocate(coord(3,nat,istep))
allocate(coord_sh(3,nat,istep))
allocate(COM_x(istep))
allocate(COM_y(istep))
allocate(COM_z(istep))

! Read in coordinates at each time step from OUTCAR
k=istep

do while (k.gt.0)
 read (7,'(1X,A69)') t
! Search for each new configuration for sequential time steps
 if (t=='POSITION                                       TOTAL-FORCE (eV/Angst)') then
   read(7,*)
   do i=1,nat
    read(7,*) coord(1,i,istep-k+1),coord(2,i,istep-k+1), coord(3,i,istep-k+1)
   enddo
   k=k-1
 endif
enddo

close(7)

! Calculate center of mass of configuration for each time step
do k=1,istep
   COM_x(k) = 0.0
   COM_y(k) = 0.0
   COM_z(k) = 0.0
   do i=1,natom1
      COM_x(k) = COM_x(k) + mass1*coord(1,i,k)
      COM_y(k) = COM_y(k) + mass1*coord(2,i,k)
      COM_z(k) = COM_z(k) + mass1*coord(3,i,k)
   enddo
   do i=natom1+1,natom1+natom2
      COM_x(k) = COM_x(k) + mass2*coord(1,i,k)
      COM_y(k) = COM_y(k) + mass2*coord(2,i,k)
      COM_z(k) = COM_z(k) + mass2*coord(3,i,k)
   enddo
   do i=natom1+natom2+1,natom1+natom2+natom3
      COM_x(k) = COM_x(k) + mass3*coord(1,i,k)
      COM_y(k) = COM_y(k) + mass3*coord(2,i,k)
      COM_z(k) = COM_z(k) + mass3*coord(3,i,k)
   enddo
COM_x(k) = COM_x(k)/mass_tot
COM_y(k) = COM_y(k)/mass_tot
COM_z(k) = COM_z(k)/mass_tot
enddo

open(9,file='COM',status='new')

do k=1,istep
   write(9,*) COM_x(k), COM_y(k), COM_z(k)
enddo
close(9)

! Create file for each configuration for each time step
do i =1, istep
 if (i<10) then
  tmp='(A3,I1)'
 elseif ((i.ge.10).and.(i.lt.100)) then
  tmp='(A3,I2)'
 elseif ((i.ge.100).and.(i.lt.1000)) then
  tmp='(A3,I3)'
 else
  tmp='(A3,I4)'
 endif

! Write out COM-shifted coordinates for each time step X to posX
 write(filename,tmp) 'pos',i     ! filename=pos+char(i)

 open(8,file=filename,status='new')
 rewind(8)
 do k=1,natom1
    coord_sh(1,k,i) = coord(1,k,i)-COM_x(i)
    coord_sh(2,k,i) = coord(2,k,i)-COM_y(i)
    coord_sh(3,k,i) = coord(3,k,i)-COM_z(i)
    write(8,'(1X,A2,1X,6(F12.5,1X))') atom1,coord_sh(1,k,i),coord_sh(2,k,i),coord_sh(3,k,i)
 enddo
 do k=natom1+1,natom1+natom2
    coord_sh(1,k,i) = coord(1,k,i)-COM_x(i)
    coord_sh(2,k,i) = coord(2,k,i)-COM_y(i)
    coord_sh(3,k,i) = coord(3,k,i)-COM_z(i)
    write(8,'(1X,A2,1X,6(F12.5,1X))') atom2,coord_sh(1,k,i),coord_sh(2,k,i),coord_sh(3,k,i)
 enddo
 do k=natom1+natom2+1,natom1+natom2+natom3
    coord_sh(1,k,i) = coord(1,k,i)-COM_x(i)
    coord_sh(2,k,i) = coord(2,k,i)-COM_y(i)
    coord_sh(3,k,i) = coord(3,k,i)-COM_z(i)
    write(8,'(1X,A2,1X,6(F12.5,1X))') atom3,coord_sh(1,k,i),coord_sh(2,k,i),coord_sh(3,k,i)
 enddo
 close(8)
enddo


stop
end program

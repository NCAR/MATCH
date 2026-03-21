c     $Header$

c     Purpose: library of Fortran string manipulation routines

c     Usage: 
c     cd /home/zender/csz; make -W sng.F OPTS=D libcsz; cd -
c     cd /home/zender/csz; make -W sng.F libcsz; cd -
c     cd /home/zender/csz; make OPTS=D libcsz; cd -

c     NB: All Fortran compilers allow dimensioning of subroutine arguments by other subroutine arugments
c     In addition, g77 allows dimensioning of local variables by subroutine arguments (i.e., wvl_out_nbr)
c     Sun f77, however, requires that local variables be dimensioned by fixed parameters

      integer function f_strlen(sng)
c     Purpose: Return the length of the string preceding the first null character
c     Fortran intrinsic len(sng) is returned if the string is not null terminated
c     Prototype:
c     integer f_strlen            ! [nbr] Length of string (libcsz_f77)
      implicit none
c     Input
      character*(*) sng
c     Local
      integer idx
      integer lng

c     Main Code
      lng=len(sng)
      do idx=1,lng
         if (sng(idx:idx) == char(0)) goto 100
      end do                    ! end loop over characters
 100  f_strlen=idx-1
      return
      end                       ! end f_strlen()

      integer function f_strstr(sng1,sng2)
c     Purpose: Return the location of second string in the first
c     Fortran intrinsic len(sng) is returned if the string is not null terminated
c     Prototype:
c     integer f_strstr            ! Location of sng1 in sng2 (libcsz_f77)
      implicit none
c     Input
      character*(*) sng1
      character*(*) sng2
c     Local
      integer len1
      integer len2
      integer lng1
      integer lng2
      integer idx1
      integer idx2
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)

c     Main Code
      len1=len(sng1)
      len2=len(sng2)
      lng1=f_strlen(sng1)
      lng2=f_strlen(sng2)
      do idx1=1,lng1
         if (sng1(idx1:idx1) == sng2(1:1)) then
            do idx2=2,lng2
               if (sng1(idx1+idx2-1:idx1+idx2-1) /= sng2(idx2:idx2)) goto 200
            end do              ! end loop over characters
            goto 100
         endif                  ! endif
 200     continue
      end do                    ! end loop over characters
 100  continue
      if (idx2 == lng2+1) then
         f_strstr=idx1
      else
         f_strstr=-1
      endif                     ! endif
      return
      end                       ! end f_strstr()

      subroutine strini(sng)
c     Purpose: Initialize all elements of a character array to char(0)
      implicit none
c     Input/Output
      character*(*) sng
c     Local
      integer idx
      integer lng

c     Main Code
      lng=len(sng)
      do idx=1,lng
         sng(idx:idx)=char(0)
      end do                    ! end loop over characters
c      write (6,*) 'strini(): Initialized string of length',lng
      return
      end                       ! end strini()

      subroutine strnul(sng)
c     Purpose: Change space characters to null characters in a string
c     Only spaces after the last significant character in string are affected.
c     Spaces before last significant character in string are not affected.
      implicit none
c     Input/Output
      character*(*) sng
c     Local
      integer idx
      integer lng
      integer lst_sgn_chr_idx
c     External
      integer strlsc            ! Last significant character (libcsz_f77)
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)

c     Main Code
      lng=f_strlen(sng)
      lst_sgn_chr_idx=strlsc(sng)
      do idx=lst_sgn_chr_idx+1,lng
         if(sng(idx:idx) == ' ') sng(idx:idx)=char(0)
      end do                    ! end loop over characters
      return
      end                       ! end strnul()

      character*(*) function strnulc(sng)
c     Purpose: Change space characters to null characters in a string
c     Only spaces after the last significant character in string are affected.
c     Spaces before last significant character in string are not affected.
      implicit none
c     Input/Output
      character*(*) sng
c     Local
      integer idx
      integer lng
      integer lst_sgn_chr_idx
c     External
      integer strlsc            ! Last significant character (libcsz_f77)
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)

c     Main Code
      lng=f_strlen(sng)
      lst_sgn_chr_idx=strlsc(sng)
      do idx=lst_sgn_chr_idx+1,lng
         if(sng(idx:idx) == ' ') sng(idx:idx)=char(0)
      end do                    ! end loop over characters
      strnulc(1:lng)=sng(1:lng)
      return
      end                       ! end strnulc()

      integer function strlsc(sng)
c     Purpose: Return position of last significant character in string
      implicit none
c     Input
      character*(*) sng
c     Local
      integer idx
      integer lng

c     Main Code
      lng=len(sng)
      do idx=lng,1,-1
c     Look for last normal character
         if (sng(idx:idx) /= ' '.and.sng(idx:idx) /= char(0)) goto 100
      end do                    ! end loop over characters
 100  continue
      strlsc=idx
      return
      end                       ! end strlsc()

      subroutine f_strcpy(sng1,sng2)
c     Purpose: Copy sng2 into sng1
c     Any space remaining at the end of sng1 is null-initialized
c     This function works well to at copying fixed strings into variables, e.g.,
c     call f_strcpy(sng,'Hello World')
      implicit none
c     Commons
#include <dbg.com> /* Debugging constants */
c     Input
      character*(*) sng2
c     Input/Output
      character*(*) sng1
c     Local
      integer len1
      integer len2
      integer lng1
      integer lng2
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)

c     Main Code
      len1=len(sng1)
      len2=len(sng2)
      lng1=f_strlen(sng1)
      lng2=f_strlen(sng2)
      if (len1 < lng2) then
         write (6,'(a,a)') prg_nm(1:f_strlen(prg_nm)),': ERROR len1 < lng2 in f_strcpy()'
         write (6,'(a,i3,a,i3,2a)') 'len1 = ',len1,', lng1 = ',lng1,', sng1 = ',sng1
         write (6,'(a,i3,a,i3,2a)') 'len2 = ',len2,', lng2 = ',lng2,', sng2 = ',sng2
         stop 'EXIT_FAILURE from f_strcpy()'
      endif                     ! endif
      sng1(1:lng2)=sng2(1:lng2)
c     Null-initialize rest of string
      call strnul(sng1)
      return
      end                       ! end f_strcpy()

      subroutine f_strcat(sng1,sng2)
c     Purpose: sng1 := sng1 // sng2
c     The case where sng1=sng2 is handled correctly
      implicit none
c     Commons
#include <dbg.com> /* Debugging constants */
c     Input
      character*(*) sng2
c     Output
      character*(*) sng1
c     Local
      integer len1
      integer len2
      integer lng1
      integer lng2
c     External
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)

c     Main Code
      len1=len(sng1)
      len2=len(sng2)
      lng1=f_strlen(sng1)
      lng2=f_strlen(sng2)
      if (lng1+lng2 >= len1) then
         write (6,'(a,a)') prg_nm(1:f_strlen(prg_nm)),': ERROR lng1+lng2 >= len1 in f_strcat()'
         write (6,'(a,i3,a,i3,2a)') 'len1 = ',len1,', lng1 = ',lng1,', sng1 = ',sng1
         write (6,'(a,i3,a,i3,2a)') 'len2 = ',len2,', lng2 = ',lng2,', sng2 = ',sng2
         stop 'EXIT_FAILURE from f_strcat()'
      endif                     ! endif
      sng1(lng1+1:lng1+lng2)=sng2(1:lng2)
      return
      end                       ! end f_strcat()

      subroutine prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID)
c     Purpose: Create identity string from CVS input
      implicit none
c     Commons
#include <dbg.com> /* Debugging constants */
c     Parameter
      integer CVS_kk
      integer CVS_kv
      integer CVS_kkv
      parameter(CVS_kk=0)
      parameter(CVS_kv=1)
      parameter(CVS_kkv=2)
c     Input
      character*(*) CVS_Id
      character*(*) CVS_Revision
      character*(*) CVS_Date
c     Output
      character*(*) prg_ID
c     Local
      integer CVS_typ
      integer Date_ptr
      integer slash_ptr
      integer prg_ptr
      integer vrs_ptr
c     External
      integer f_strstr
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)
c     Main Code
c     Decide how CVS keywords were expanded
      Date_ptr=f_strstr(CVS_Date,'Date')
      if (Date_ptr < 0) then
c     String 'Date' was not found---CVS expansion must be -kv
         CVS_typ=CVS_kv
      else
c     String 'Date' was found---CVS expansion is -kk or -kkv
         slash_ptr=f_strstr(CVS_Date,'/')
         if (slash_ptr > 0) then 
            CVS_typ=CVS_kkv 
         else 
            CVS_typ=CVS_kk
         endif                  ! endif
      endif                     ! endif

      if (.false.) then
         write (6,'(a,a)') prg_nm(1:f_strlen(prg_nm)),': DEBUG prg_ID_mk() reports CVS strings:'
         write (6,'(a,a)') 'CVS_Id: ',CVS_Id
         write (6,'(a,a)') 'CVS_Revision: ',CVS_Revision
         write (6,'(a,a)') 'CVS_Date: ',CVS_Date
         write (6,'(a,i1)') 'CVS_typ: ',CVS_typ
         write (6,'(a,i2)') 'Date_ptr: ',Date_ptr
         write (6,'(a,i2)') 'slash_ptr: ',slash_ptr
      endif

      call strnul(CVS_Revision)
      call strnul(CVS_Date)
      call strini(prg_ID)
      if (CVS_typ == CVS_kk) then
         prg_ID='Source file unknown Version unknown Date unknown' // char(0)
      else 
         prg_ptr=f_strstr(CVS_Id,',v') ! ',v' is right after program name
         if (prg_ptr < 0) stop 'ERROR: prg_ID_mk() unable to find source code name'
         if (CVS_typ == CVS_kv) then
            prg_ID=CVS_Id(1:prg_ptr-1) // 
     $           ' version ' // CVS_Revision(1:f_strlen(CVS_Revision)) // 
     $           ' dated ' // CVS_Date(1:f_strlen(CVS_Date)) // ' GMT' // char(0)
         else if (CVS_typ == CVS_kkv) then
            vrs_ptr=f_strstr(CVS_Revision,'ion:') ! 'ion:' should appear two characters before version
            if (vrs_ptr < 0) stop 'ERROR: prg_ID_mk() unable to find revision number'
            prg_ID=CVS_Id(6:prg_ptr-1) // 
     $           ' version ' // CVS_Revision(vrs_ptr+5:f_strlen(CVS_Revision)-2) // 
     $           ' dated ' // CVS_Date(8:26) // ' GMT' // char(0)
         else 
            write (6,'(a,a)') prg_nm(1:f_strlen(prg_nm)),': ERROR prg_ID_mk() reports unknown CVS_typ'
         endif                  ! endif
      endif                     ! endif
      return
      end                       ! end prg_ID_mk()
      
      subroutine cmd_ln_sng(cmd_ln)
c     Purpose: Return a copy of the command line and initialize the program name
      implicit none
c     Commons
#include <dbg.com> /* Debugging constants */
c     Output
      character cmd_ln*(*)  
c     Local
      character argv*80         ! Maximum length of any single argument
      integer arg_idx
      integer arg_nbr
      integer cmd_ln_len
      integer cmd_ln_lng
c     External
      integer iargc
      integer f_strlen            ! [nbr] Length of string (libcsz_f77)

c     Main Code

c     Initialize defaults
      call strini(cmd_ln)
      call strini(prg_nm)
      cmd_ln_len=len(cmd_ln)
      
      arg_nbr=iargc()
      if (dbg_lvl >= dbg_vec) then
         write (6,'(a,i2)') 'cmd_ln_sng(): arg_nbr = ',arg_nbr
      endif                     ! endif dbg

c     Loop over arguments
      do arg_idx=0,arg_nbr
         if (arg_idx > 0) call f_strcat(cmd_ln,' ')
         call strini(argv)
         call getarg(arg_idx,argv)
         call strnul(argv)
         call f_strcat(cmd_ln,argv)
         if (arg_idx == 0) call f_strcpy(prg_nm,argv)
         cmd_ln_lng=f_strlen(cmd_ln)
         if (cmd_ln_lng > cmd_ln_len) stop 'cmd_ln_lng > cmd_ln_len in cmd_ln_sng()'
         if (dbg_lvl >= dbg_vrb) then
            write (6,'(a,i3)') 'arg_idx = ',arg_idx
            write (6,'(a,a)') 'argv = ',argv
            write (6,'(a,i3)') 'len(argv) = ',len(argv)
            write (6,'(a,i3)') 'f_strlen(argv) = ',f_strlen(argv)
            write (6,'(a,i3)') 'f_strlen(cmd_ln) = ',f_strlen(cmd_ln)
            write (6,'(a,a)') 'cmd_ln = ',cmd_ln
         endif                  ! endif dbg
      enddo                     ! end loop over arg

      if (dbg_lvl >= dbg_io) write (6,'(a,a)') 'cmd_ln = ',cmd_ln(1:f_strlen(cmd_ln))
      return
      end                       ! end cmd_ln_sng()

      character*10 function date2sng(idate)
c     Purpose: Convert integer date in YYYYMMDD to a character array in the form
c     'YYYY-MM-DD'
c     Author: Brian Eaton charutl.F:idate2char()
c     Recoded: 19990901 Charlie Zender
c     Usage:
c     character*10 date2sng   ! [sng] Convert YYYYMMDD integer to YYYY-MM-DD string (libcsz_f77)
      implicit none
c     Input
      integer idate
c     Local
      integer yy 
      integer mm
      integer dd
      character dash
c     Main code
      dash='-'
c     Extract year, month, and day from date
      yy=idate/10000
      mm=mod(abs(idate),10000)/100
      dd=mod(abs(idate),100)
      write(date2sng,'(i4.4,a1,i2.2,a1,i2.2)') yy,dash,mm,dash,dd
      return
      end

      character*8 function sec2sng(isec)
c     Purpose: Convert integer seconds to a character array in the form 'HH:MM:SS'
c     Author: Brian Eaton charutl.F:isec2char()
c     Recoded: 19990901 Charlie Zender
c     character*8 sec2sng   ! [sng] Convert integer to 'HH:MM:SS' string (libcsz_f77)
      implicit none
c     Input
      integer isec
c     Local
      integer hh
      integer mm
      integer ss
      character colon
c     Main Code
      colon=':'
c     Extract hour, minutes, seconds
      hh=isec/3600
      mm=mod(isec,3600)/60
      ss=mod(isec,60)
      write(sec2sng,'(i2.2,a1,i2.2,a1,i2.2)') hh,colon,mm,colon,ss
      return
      end

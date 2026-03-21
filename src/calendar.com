c $Id$

      character(len=16) ::
     $  type  ! calendar type.  Currently '365' or 'gregorian'

      common /calendar_C/ type

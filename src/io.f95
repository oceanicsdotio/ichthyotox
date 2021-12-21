module io

    use variables, only : ZERO, sp
    use random, only : random_number_generator
    
    implicit none

  contains
  
    subroutine getInteger(file, key, target, fid)
  
      character(len = *), intent(in) :: key, file
      integer, intent(out) :: target
      integer, intent(in), optional :: fid
      integer :: errorCode
  
      errorCode = find_key(file, key, iscal=target)
      if (errorCode /= 0) then
        write(fid, *) 'Error reading '//key//': ', errorCode
        stop
      end if
  
    end subroutine
  
  
    subroutine getString(file, key, target, fid)
  
      character(len = *), intent(in) :: key, file
      character(len = 80), intent(out) :: target
      integer, intent(in) :: fid
      integer :: errorCode, ii
  
      errorCode = find_key(file, key, cval=target)
      if (errorCode /= 0) then
        write(fid, *) 'Error reading '//key//': ', errorCode
        stop
      else
        ! remove trailing directory '/'
        ii = len_trim(target)
        if (target(ii:ii) == "/") target(ii:ii) = " "
      end if
  
    end subroutine
  
  
    integer function find_key(FNAME, VNAME, ISCAL, FSCAL, IVEC, FVEC, CVEC, NSZE, CVAL, LVAL)
      !   Scan an Input File for a Variable
      !   RETURN VALUE:
      !        0 = FILE FOUND, VARIABLE VALUE FOUND
      !       -1 = FILE DOES NOT EXIST OR PERMISSIONS ARE INCORRECT
      !       -2 = VARIABLE NOT FOUND OR IMPROPERLY SET
      !       -3 = VARIABLE IS OF DIFFERENT TYPE, CHECK INPUT FILE
      !       -4 = VECTOR PROVIDED BUT DATA IS SCALAR TYPE
      !       -5 = NO DATATYPE DESIRED, EXITING
  
      !   REQUIRED INPUT:
      !        FNAME = File Name
      !        FSIZE = Length of Filename
  
      !   optional (MUST PROVIDE ONE)
      !        ISCAL = integer SCALAR
      !        FSCAL = FLOAT SCALAR
      !        CVAL = character VARIABLE
      !        LVAL = LOGICAL VARIABLE
      !        IVEC = integer VECTOR **
      !        FVEC = FLOAT VECTOR **
      !        CVEC = STRING VECTOR **
      !      **NSZE = ARRAY SIZE (MUST BE PROVIDED WITH IVEC/FVEC)
  
      use variables
  
      character(LEN = *) :: FNAME, VNAME
      integer, intent(inout), optional :: ISCAL, IVEC(*)
      REAL(SP), intent(inout), optional :: FSCAL, FVEC(*)
      character(LEN=80), optional :: CVAL, CVEC(*)
      LOGICAL, intent(inout), optional :: LVAL
      integer, intent(inout), optional :: NSZE
  
      REAL(SP) REALVAL(150)
      integer  INTVAL(150)
      character(LEN=20 ) :: key
      character(LEN=80 ) :: STRINGVAL(150),TITLE
      character(LEN=80 ) :: line
      character(LEN=400) :: buffer
      character(LEN=7  ) :: type_of
      character(LEN=20 ), DIMENSION(200) :: SET
      integer :: last, NVAL, lines, NREP
      logical :: SETYES, ALLSET, CHECK, LOGVAL
      character(len=*), parameter :: continue_line = "////"
      character(len = len_trim(copy)) :: text
      character(len = len_trim(copy)) :: value, TEMP, fragments(200)
      character(len = 80) :: TSTRING
      character(len = 6) :: ERRSTRING
      character(len = 16) :: NUMCHARS = "0123456789+-Ee. "
      integer :: EQLOC, length, ii, LOCEX, NP
      logical :: flag
  
      find_key = 0
  
      ! OPEN THE INPUT FILE
      inquire(file=TRIM(FNAME), exist=CHECK)
      if (.not. CHECK) then
        find_key = -1
      end if
  
      open(10, file=trim(FNAME))
      rewind(10)
  
      lines = 0
      do while (.true.)
  
        buffer(1:len(buffer)) = ' '
        NREP  = 0
        lines = lines + 1
        read(10,'(a)', end=20) line
        buffer(1:80) = line(1:80)
  
        ! PROCESS LINE CONTINUATIONS
        110 CONTINUE
        last = len_trim(line)
        if (last /= 0) then
          if ( line(last-1:last) == '\\\\') then
  
            NREP = NREP + 1
            read(10, '(a)', end=20) line
            lines = lines + 1
            buffer( NREP*80 + 1 : NREP*80 +80) = line(1:80)
            GOTO 110
          end if
        end if
  
        ! REMOVE LINE CONTINUATION character \\
        if (NREP > 0) then
          do last = 2, LEN_TRIM(buffer)
            if ( buffer(last-1:last) == '\\\\') then
              buffer(last-1:last) = '  ' 
            endif
          enddo
        endif
  
        fragments = " "
        type_of = "error"
        LOGVAL = .false.
  
        write(ERRSTRING, "(I6)") lines
  
        LOCEX = index(text, "!")
        if (LOCEX /= 0) text = text(1:LOCEX-1)
        length = len_trim(text)
  
        if (length == 0) then
          type_of = "none"
          key = "none"
          return
        end if
  
        ! Commas to spaces
        where (text == ",")
          text = " "
        end where
  
        ! Find assignment "="
        EQLOC = index(text, "=")
        if (EQLOC == 0) call raise(6,'DATA LINE '//ERRSTRING//' MUST CONTAIN "=" ')
  
        ! split name and value substrings
        key = text(1:EQLOC-1)
        value  = adjustl(text(EQLOC+1:LENGTH))
        length = len_trim(value)
  
        if (length == 0) call raise(6,'IN DATA PARAMETER FILE', 'VARIABLE LINE'//ERRSTRING//' HAS NO ASSOCIATED VALUE')
  
        ! check for logical
        if ((value(1:1) == "T" .or. value(1:1) == "F") .and. length == 1) then
          type_of = "logical"
          if (value(1:1) == "T") LOGVAL = .true.
          return
        end if
  
        ! is string if contains non-numeric characters
        do ii = 1, length
          if (index(NUMCHARS, value(ii:ii)) == 0) then
  
            type_of = "string"
            TSTRING = value
            stringval(1) = TSTRING
            NVAL = 1
            flag = .true.
  
            do ii = 1, length
              if (value(ii:ii) /= " ") then
                fragments(NVAL) = trim(fragments(NVAL)) // value(ii:ii)
                flag = .true.
              else
                if (flag) NVAL = NVAL + 1
                flag = .false.
              end if
            end do
  
            do ii = 1, NVAL
              stringval(ii + 1) = trim(fragments(ii))
            end do
            return
  
          end if
        end do
  
        type_of = merge("float  ", "integer", index(value, ".") /= 0)
  
        ! Split lines
        NP = 1
        flag = .true.
        do ii = 1, length
          if (value(ii:ii) /= " ") then
            fragments(NP) = trim(fragments(NP)) // value(ii:ii)
            flag = .true.
          else
            if (flag) NP = NP + 1
            flag = .false.
          end if
        end do
  
        ! numerical
        NVAL = NP
        do ii = 1, NP
          if (type_of == "float") then
            read(trim(fragments(ii)), *) realval(ii)
          else
            read(trim(fragments(ii)), *) intval(ii)
          end if
        end do
  
  
        if (trim(name) == trim(VNAME)) then
  
          if (PRESENT(ISCAL)) then
            if (type_of == 'integer') then
              ISCAL = INTVAL(1)
              return
            else
              find_key = -3
            end if
          elseif(present(FSCAL)) then
            if (type_of == 'float') then
              FSCAL = REALVAL(1)
              return
            else
              find_key = -3
            end if
          elseif(present(CVAL))THEN
            if (type_of == 'string') then
              CVAL = STRINGVAL(1)
              return
            else
              find_key = -3
            end if
          elseif (present(LVAL)) THEN
            if (type_of == 'logical') then
              LVAL = LOGVAL
              return
            else
              find_key = -3
            end if
          else if (present(IVEC)) then
            if (NVAL > 1) then
              if (type_of == 'integer') then
                IVEC(1:NVAL) = INTVAL(1:NVAL)
                NSZE = NVAL
                return
              else
                find_key = -3
              end if
            else
              find_key = -4
            end if
          elseif (present(FVEC)) then
            if (NVAL > 1) then
              IF (type_of == 'float') then
                FVEC(1:NVAL) = REALVAL(1:NVAL)
                NSZE = NVAL
                return
              else
                find_key = -3
              end if
            else
              find_key = -4
            end if
          elseif (present(CVEC)) then
            if (NVAL > 0) then
              if (type_of == 'string') then
                CVEC(1:NVAL) = STRINGVAL(2:NVAL+1)
                NSZE = NVAL
                return
              else
                find_key = -3
              end if
            else
              find_key = -4
            end if
          else
            find_key = -5
          end if
        end if
      end do
      20 close(10)
      find_key = -2
    end function
  
  end module
  
module MOD_INP

contains

    subroutine DATA_RUN
        ! Input Parameters Which Control the Model Run
        use ALL_VARS

        implicit none

        integer :: ISCAN
        character(len = 120) :: filename
        integer :: ii


        ! Read in variables and set values
        filename = "./"//trim(folderprefix)//"/"//trim(CASENAME)//"_run.dat"

        ! Info file
        ISCAN = scan(filename, "INFOFILE", CVAL = INFOFILE)
        if (ISCAN /= 0) then
            write(IPT, *) 'ERROR READING INFOFILE: ', ISCAN
            stop
        end if

        ! Open runtime info file
        IPT = 71
        if (trim(INFOFILE) /= "screen") then
            open(IPT, FILE = trim(INFOFILE))
        else
            IPT = 6
        end if

        ! External time step (DTI)
        ISCAN = scan(trim(filename), "DTI", FSCAL = DTI)
        if (ISCAN /= 0) then
            write (IPT, *) 'ERROR READING DTI: ', ISCAN
            stop
        end if

        ! Input time step of flow fields (instp)
        ISCAN = scan(trim(filename),"INSTP", FSCAL = INSTP)
        if (ISCAN /= 0) then
            write(IPT, *) 'ERROR READING INSTP: ', ISCAN
            stop
        end if

        ! External time step (DTOUT)
        ISCAN = scan(trim(filename), "DTOUT", FSCAL = DTOUT)
        if (ISCAN /= 0) then
            write (IPT, *) 'ERROR READING DTOUT: ', ISCAN
            stop
        end if

        ! Total time to move drifters (TDRIFT)
        ISCAN = scan(trim(filename), "TDRIFT", ISCAL = TDRIFT)
        if (ISCAN /= 0)then
            write(IPT, *) 'ERROR READING TDRIFT: ', ISCAN
            stop
        end if


        ! Input year of run (YEARLAG)
        ISCAN = scan(trim( filename),"YEARLAG",ISCAL = YEARLAG)
        if (ISCAN /= 0) then
            write(IPT, *)'ERROR READING YEARLAG: ', ISCAN
            stop
        end if

        ! Input month of run (MONTHLAG)
        ISCAN = scan(trim(filename), "MONTHLAG", ISCAL = MONTHLAG)
        if (ISCAN /= 0) then
            write(IPT, *) 'ERROR READING YEARLAG: ', ISCAN
            stop
        end if

        ! Input day of run (DAYLAG)
        ISCAN = scan(trim(filename), "DAYLAG", ISCAL = DAYLAG)
        if (ISCAN /= 0) then
            write(IPT, *) 'ERROR READING DAYLAG: ', ISCAN
            stop
        end if

        ! Input hour of run (HOURLAG)
        ISCAN = scan(trim(filename), "HOURLAG", ISCAL = HOURLAG)
        if (ISCAN /= 0) then
            write(IPT, *) 'ERROR READING HOURLAG: ', ISCAN
            stop
        end if

        ! "P_SIGMA" turns on vertical location of particles in sigma
        ISCAN = scan(trim(filename), "P_SIGMA", LVAL = P_SIGMA)
        if (ISCAN /= 0) then
            write(IPT, *) 'ERROR READING P_SIGMA: ', ISCAN
            if (ISCAN == -2) then
                write(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',trim(filename)
            end if
            stop
        end if

        ! "OUT_SIGMA" TURNS ON VERTICAL LOCATION OF PARTICLES IN SIGMA
        ISCAN = scan(trim(filename), "OUT_SIGMA", LVAL = OUT_SIGMA)
        if (ISCAN /= 0) then
            write(IPT, *) 'ERROR READING OUT_SIGMA: ', ISCAN
            if (ISCAN == -2) then
                write(IPT, *) 'VARIABLE NOT FOUND IN INPUT FILE: ', trim(filename)
            end if
            stop
        end if


        ! "F_DEPTH" KEEP SAME Z DEPTH ALONG THE TRACKING
        ISCAN = scan(trim(filename), "F_DEPTH", LVAL = F_DEPTH)
        if (ISCAN /= 0) then
            write(IPT, *) 'ERROR READING F_DEPTH: ', ISCAN
            if (ISCAN == -2) then
                write(IPT, *) 'VARIABLE NOT FOUND IN INPUT FILE: ', trim(filename)
            end if
            stop
        end if

        ! RANDOM WALK CHOICE
        ISCAN = scan(trim(filename), "IRW", ISCAL = IRW)
        if (ISCAN /= 0) then
            write(IPT, *) 'ERROR READING IRW: ', ISCAN
            stop
        end if

        ! Horizontal diffusion coefficient (DHOR)
        ISCAN = scan(trim(filename), "DHOR", FSCAL = DHOR)
        if (ISCAN /= 0) then
            write(IPT, *) 'ERROR READING DHOR: ', ISCAN
            stop
        end if

        !  RANDOM WALK TIME STEP (DTRW)
        ISCAN = scan(trim(filename), "DTRW", FSCAL = DTRW)
        if (ISCAN /= 0) then
            write(IPT, *) 'ERROR READING DTRW: ', ISCAN
            stop
        end if

        ! "GEOAREA" DIRECTORY FOR INPUT FILES
        ISCAN = scan(filename,"GEOAREA",CVAL = GEOAREA)
        if (ISCAN /= 0) then
            write(IPT, *) 'ERROR READING GEOAREA: ', ISCAN
            stop
        end if

        ii = len_trim(GEOAREA)
        if(GEOAREA(ii:ii) == "/") GEOAREA(ii:ii) = " "


        ! DIRECTORY FOR INPUT FILES
        ISCAN = scan(filename, "INPDIR", CVAL = INPDIR)
        if (ISCAN /= 0) then
            write(IPT, *)'ERROR READING INPDIR: ', ISCAN
            stop
        end if

        ii = len_trim(INPDIR)
        if(INPDIR(ii:ii) == "/") INPDIR(ii:ii) = " "


        ! INPUT FILES
        ISCAN = scan(filename, "LAGINI", CVAL = LAGINI)
        if (ISCAN /= 0) then
            write(IPT, *) 'ERROR READING LAGINI: ', ISCAN
            stop
        end if

        ii = len_trim(LAGINI)
        if (LAGINI(ii:ii) == "/") LAGINI(ii:ii) = " "


        ! "OUTDIR"
        ISCAN = scan(filename, "OUTDIR", CVAL = OUTDIR)
        if (ISCAN /= 0) then
            write(IPT, *) 'ERROR READING OUTDIR: ', ISCAN
            stop
        end if

        ii = len_trim(OUTDIR)
        if (OUTDIR(ii:ii) == "/") OUTDIR(ii:ii) = " "

        ! Set unit values for input output files
        IOPAR=11
        INLAG=13

    end subroutine

    subroutine parse_line(LNUM, NUMCHAR, key, type_of, LOGVAL, STRINGVAL, REALVAL, INTVAL, NVAL)
        ! Decompose input line into variable name and value
        use parameters
        implicit none

        intrinsic adjustl

        integer, intent(in) :: LNUM, NUMCHAR
        character(len = NUMCHAR) :: text
        character(len = 20), intent(out) :: key
        character(len = 7), intent(out) :: type_of

        logical, intent(out) :: LOGVAL
        character(LEN = 80), intent(out), dimension(150) :: STRINGVAL
        real(SP), intent(inout), dimension(150) :: REALVAL
        integer, intent(inout), dimension(150) :: INTVAL
        integer, intent(out) :: NVAL

        character(len = NUMCHAR) :: value, TEMP, fragments(200)
        character(len = 80) :: TSTRING
        character(len = 6) :: ERRSTRING
        character(len = 16) :: NUMCHARS
        integer :: EQLOC, length, ii, LOCEX, NP
        logical :: flag

        fragments = " "
        NUMCHARS = "0123456789+-Ee. "
        type_of = "error"
        LOGVAL = .false.
        
        write(ERRSTRING, "(I6)") LNUM
        
        LOCEX = index(text, "!")
        if (LOCEX /= 0) text = text(1:LOCEX-1)
        length = len_trim(text)
        
        if (length == 0) then
            type_of = "none"
            key = "none"
            return
        end if

        ! Commas to spaces
        where (text(:) == ",")
            text(:) = " "
        end where
        
        ! Find assignment "="
        EQLOC = INDEX(text,"=")
        IF (EQLOC == 0) CALL raise(6,'DATA LINE '//ERRSTRING//' MUST CONTAIN "=" ')

        ! split name and value substrings
        key = text(1:EQLOC-1)
        value  = adjustl(text(EQLOC+1:LENGTH))
        length = len_trim(value)
        
        IF (length == 0) CALL raise(6,'IN DATA PARAMETER FILE', 'VARIABLE LINE'//ERRSTRING//' HAS NO ASSOCIATED VALUE')

        ! DETERMINE TYPE OF value
        ! CHECK FOR LOGICAL
        IF (((value(1:1) == "T") .or. (value(1:1) == "F")) .and. (length == 1)) then
            type_of = "logical"
            if (value(1:1) == "T") LOGVAL = .true.
            return
        end if

        ! CHECK IF IT IS A STRING  (CONTAINS characterS OTHER THAN 0-9,+,-,e,E,.)
        do ii = 1, length
            if (INDEX(NUMCHARS, value(ii:ii)) == 0) type_of = "string"
        end do

        ! PROCESS STRING (MAY BE MULTIPLE)
        if (type_of == "string") then
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

        type_of = merge("float", "integer", index(value, ".") /= 0)
        
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
            read(trim(fragments(ii)), *) merge(realval(ii), intval(ii), type_of == "float")
        end do

    end subroutine
    

    integer function scan(FNAME, VNAME, ISCAL, FSCAL, IVEC, FVEC, CVEC, NSZE, CVAL, LVAL)
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

        use parameters
        implicit none
        character(LEN = *) :: FNAME, VNAME
        integer, intent(inout), optional :: ISCAL, IVEC(*)
        REAL(SP), intent(inout), optional :: FSCAL, FVEC(*)
        character(LEN=80), optional :: CVAL, CVEC(*)
        LOGICAL, intent(inout), optional :: LVAL
        integer, intent(inout), optional :: NSZE

        REAL(SP) REALVAL(150)
        integer  INTVAL(150)
        character(LEN=20 ) :: name
        character(LEN=80 ) :: STRINGVAL(150),TITLE
        character(LEN=80 ) :: line
        character(LEN=400) :: copy
        character(LEN=7  ) :: type_of
        character(LEN=20 ), DIMENSION(200) :: SET
        integer :: last, NVAL, lines, NREP
        logical :: SETYES, ALLSET, CHECK, LOGVAL
        character(len=*), parameter :: continue_line = "////"


        scan = 0

        ! OPEN THE INPUT FILE
        inquire(file=TRIM(FNAME), exist=CHECK)
        if (.not. CHECK) then
            scan = -1
        end if

        open(10, file=trim(FNAME))
        rewind(10)
        
        lines = 0
        do while (.true.)

            copy(1:len(copy)) = ' '
            NREP  = 0
            lines = lines + 1
            read(10,'(a)', end=20) line
            copy(1:80) = line(1:80)

            ! PROCESS LINE CONTINUATIONS
            110 CONTINUE
            last = len_trim(line)
            if (last /= 0) then
                if ( line(last-1:last) == '\\\\') then
                   
                    NREP = NREP + 1
                    read(10, '(a)', end=20) line
                    lines = lines + 1
                    copy( NREP*80 + 1 : NREP*80 +80) = line(1:80)
                    GOTO 110
                end if
            end if

            ! REMOVE LINE CONTINUATION character \\
            if (NREP > 0) then
                do last = 2, LEN_TRIM(copy)
                    if ( copy(last-1:last) == '\\\\') copy(last-1:last) = '  '
                end do
            end if

            call parse_line(lines, len_trim(copy), adjustl(copy), name, type_of, LOGVAL, STRINGVAL, REALVAL, INTVAL, NVAL)

            ! IF name MATCHES, PROCESS VARIABLE AND ERROR-CHECK

            if (trim(name) == trim(VNAME)) then

                if (PRESENT(ISCAL)) then
                    if (type_of == 'integer') then
                        ISCAL = INTVAL(1)
                        return
                    else
                        scan = -3
                    end if
                elseif(present(FSCAL)) then
                    if (type_of == 'float') then
                        FSCAL = REALVAL(1)
                        return
                    else
                        scan = -3
                    end if
                elseif(present(CVAL))THEN
                    if (type_of == 'string') then
                        CVAL = STRINGVAL(1)
                        return
                    else
                        scan = -3
                    end if
                elseif (present(LVAL)) THEN
                    if (type_of == 'logical') then
                        LVAL = LOGVAL
                        return
                    else
                        scan = -3
                    end if
                else if (present(IVEC)) then
                    if (NVAL > 1) then
                        if (type_of == 'integer') then
                            IVEC(1:NVAL) = INTVAL(1:NVAL)
                            NSZE = NVAL
                            return
                        else
                            scan = -3
                        end if
                    else
                        scan = -4
                    end if
                elseif (present(FVEC)) then
                    if (NVAL > 1) then
                        IF (type_of == 'float') then
                            FVEC(1:NVAL) = REALVAL(1:NVAL)
                            NSZE = NVAL
                            return
                        else
                            scan = -3
                        end if
                    else
                        scan = -4
                    end if
                elseif (present(CVEC)) then
                    if (NVAL > 0) then
                        if (type_of == 'string') then
                            CVEC(1:NVAL) = STRINGVAL(2:NVAL+1)
                            NSZE = NVAL
                            return
                        else
                            scan = -3
                        end if
                    else
                        scan = -4
                    end if
                else
                    scan = -5
                end if
            end if
        end do
        20 close(10)
        scan = -2
    end function

end module
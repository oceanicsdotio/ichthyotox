
subroutine DATA_RUN
    ! Input Parameters Which Control the Model Run
    use ALL_VARS
    use MOD_INP

    implicit none

    integer :: ISCAN
    character(len = 120) :: filename
    integer :: ii


    ! Read in variables and set values
    filename = "./"//trim(folderprefix)//"/"//trim(CASENAME)//"_run.dat"

    ! Info file
    ISCAN = SCAN_FILE(filename, "INFOFILE", CVAL = INFOFILE)
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
    ISCAN = SCAN_FILE(trim(filename), "DTI", FSCAL = DTI)
    if (ISCAN /= 0) then
        write (IPT, *) 'ERROR READING DTI: ', ISCAN
        stop
    end if

    ! Input time step of flow fields (instp)
    ISCAN = SCAN_FILE(trim(filename),"INSTP", FSCAL = INSTP)
    if (ISCAN /= 0) then
        write(IPT, *) 'ERROR READING INSTP: ', ISCAN
        stop
    end if

    ! External time step (DTOUT)
    ISCAN = SCAN_FILE(trim(filename), "DTOUT", FSCAL = DTOUT)
    if (ISCAN /= 0) then
        write (IPT, *) 'ERROR READING DTOUT: ', ISCAN
        stop
    end if

    ! Total time to move drifters (TDRIFT)
    ISCAN = SCAN_FILE(trim(filename), "TDRIFT", ISCAL = TDRIFT)
    if (ISCAN /= 0)then
        write(IPT, *) 'ERROR READING TDRIFT: ', ISCAN
        stop
    end if


    ! Input year of run (YEARLAG)
    ISCAN = SCAN_FILE(trim( filename),"YEARLAG",ISCAL = YEARLAG)
    if (ISCAN /= 0) then
        write(IPT, *)'ERROR READING YEARLAG: ', ISCAN
        stop
    end if

    ! Input month of run (MONTHLAG)
    ISCAN = SCAN_FILE(trim(filename), "MONTHLAG", ISCAL = MONTHLAG)
    if (ISCAN /= 0) then
        write(IPT, *) 'ERROR READING YEARLAG: ', ISCAN
        stop
    end if

    ! Input day of run (DAYLAG)
    ISCAN = SCAN_FILE(trim(filename), "DAYLAG", ISCAL = DAYLAG)
    if (ISCAN /= 0) then
        write(IPT, *) 'ERROR READING DAYLAG: ', ISCAN
        stop
    end if

    ! Input hour of run (HOURLAG)
    ISCAN = SCAN_FILE(trim(filename), "HOURLAG", ISCAL = HOURLAG)
    if (ISCAN /= 0) then
        write(IPT, *) 'ERROR READING HOURLAG: ', ISCAN
        stop
    end if

    ! "P_SIGMA" turns on vertical location of particles in sigma
    ISCAN = SCAN_FILE(trim(filename), "P_SIGMA", LVAL = P_SIGMA)
    if (ISCAN /= 0) then
        write(IPT, *) 'ERROR READING P_SIGMA: ', ISCAN
        if (ISCAN == -2) then
            write(IPT,*)'VARIABLE NOT FOUND IN INPUT FILE: ',trim(filename)
        end if
        stop
    end if

    ! "OUT_SIGMA" TURNS ON VERTICAL LOCATION OF PARTICLES IN SIGMA
    ISCAN = SCAN_FILE(trim(filename), "OUT_SIGMA", LVAL = OUT_SIGMA)
    if (ISCAN /= 0) then
        write(IPT, *) 'ERROR READING OUT_SIGMA: ', ISCAN
        if (ISCAN == -2) then
            write(IPT, *) 'VARIABLE NOT FOUND IN INPUT FILE: ', trim(filename)
        end if
        stop
    end if


    ! "F_DEPTH" KEEP SAME Z DEPTH ALONG THE TRACKING
    ISCAN = SCAN_FILE(trim(filename), "F_DEPTH", LVAL = F_DEPTH)
    if (ISCAN /= 0) then
        write(IPT, *) 'ERROR READING F_DEPTH: ', ISCAN
        if (ISCAN == -2) then
            write(IPT, *) 'VARIABLE NOT FOUND IN INPUT FILE: ', trim(filename)
        end if
        stop
    end if

    ! RANDOM WALK CHOICE
    ISCAN = SCAN_FILE(trim(filename), "IRW", ISCAL = IRW)
    if (ISCAN /= 0) then
        write(IPT, *) 'ERROR READING IRW: ', ISCAN
        stop
    end if

    ! Horizontal diffusion coefficient (DHOR)
    ISCAN = SCAN_FILE(trim(filename), "DHOR", FSCAL = DHOR)
    if (ISCAN /= 0) then
        write(IPT, *) 'ERROR READING DHOR: ', ISCAN
        stop
    end if

    !  RANDOM WALK TIME STEP (DTRW)
    ISCAN = SCAN_FILE(trim(filename), "DTRW", FSCAL = DTRW)
    if (ISCAN /= 0) then
        write(IPT, *) 'ERROR READING DTRW: ', ISCAN
        stop
    end if

    ! "GEOAREA" DIRECTORY FOR INPUT FILES
    ISCAN = SCAN_FILE(filename,"GEOAREA",CVAL = GEOAREA)
    if (ISCAN /= 0) then
        write(IPT, *) 'ERROR READING GEOAREA: ', ISCAN
        stop
    end if

    ii = len_trim(GEOAREA)
    if(GEOAREA(ii:ii) == "/") GEOAREA(ii:ii) = " "


    ! DIRECTORY FOR INPUT FILES
    ISCAN = SCAN_FILE(filename, "INPDIR", CVAL = INPDIR)
    if (ISCAN /= 0) then
        write(IPT, *)'ERROR READING INPDIR: ', ISCAN
        stop
    end if

    ii = len_trim(INPDIR)
    if(INPDIR(ii:ii) == "/") INPDIR(ii:ii) = " "


    ! INPUT FILES
    ISCAN = SCAN_FILE(filename, "LAGINI", CVAL = LAGINI)
    if (ISCAN /= 0) then
        write(IPT, *) 'ERROR READING LAGINI: ', ISCAN
        stop
    end if

    ii = len_trim(LAGINI)
    if (LAGINI(ii:ii) == "/") LAGINI(ii:ii) = " "


    ! "OUTDIR"
    ISCAN = SCAN_FILE(filename,"OUTDIR",CVAL = OUTDIR)
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
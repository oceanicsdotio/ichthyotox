MODULE MOD_INP

  CONTAINS

    SUBROUTINE GET_VAL(LNUM, NUMCHAR, TEXT_LINE, VARNAME, VARTYPE, LOGVAL, STRINGVAL, REALVAL, INTVAL, NVAL)
    ! Decompose input line into variable name and value
      USE MOD_PREC
      IMPLICIT NONE
      
      INTRINSIC ADJUSTL
      
      INTEGER, INTENT(IN) :: LNUM,NUMCHAR
      CHARACTER(LEN = NUMCHAR) :: TEXT_LINE
      CHARACTER(LEN = 20), INTENT(OUT) :: VARNAME
      CHARACTER(LEN = 7), INTENT(OUT) :: VARTYPE
      LOGICAL, INTENT(OUT) :: LOGVAL
      CHARACTER(LEN = 80), INTENT(OUT) :: STRINGVAL(150)
      REAL(SP), INTENT(INOUT) :: REALVAL(150)
      INTEGER, INTENT(INOUT) :: INTVAL(150)
      INTEGER, INTENT(OUT) :: NVAL

      CHARACTER(LEN = NUMCHAR) :: VARVAL,TEMP,FRAG(200)
      CHARACTER(LEN = 80) :: TSTRING
      CHARACTER(LEN = 6) :: ERRSTRING
      CHARACTER(LEN = 16) :: NUMCHARS 
      INTEGER :: LENGTH, EQLOC, LVARVAL, DOTLOC, ii, J, LOCEX, NP
      LOGICAL :: ONFRAG

      FRAG = " "
      NUMCHARS = "0123456789+-Ee. " 
      VARTYPE = "error"
      LOGVAL = .FALSE.
      LENGTH = LEN_TRIM(TEXT_LINE) 
      WRITE(ERRSTRING,"(I6)") LNUM
      LOCEX = INDEX(TEXT_LINE,"!")

      ! CHECK FOR BLANK LINE OR COMMENT
      IF((LENGTH .eq. 0) .OR. (LOCEX .eq. 1))THEN
        VARTYPE = "no data"
        VARNAME = "no data"
        RETURN
      END IF

      ! CHANGE COMMAS TO BLANKS
      do ii = 1, LENGTH
        if (TEXT_LINE(ii:ii) == ",") TEXT_LINE(ii:ii) = " "
      end do

    ! REMOVING TRAILING COMMENTS
      IF (LOCEX .ne. 0) THEN
        TEMP = TEXT_LINE(1:LOCEX-1)
        TEXT_LINE = TEMP
      END IF

      ! ENSURE "=" EXISTS AND DETERMINE LOCATION
      EQLOC = INDEX(TEXT_LINE,"=")
      IF(EQLOC .eq. 0) CALL PERROR(6,'DATA LINE '//ERRSTRING//' MUST CONTAIN "=" ')

      ! SPLIT OFF VARNAME AND VARVAL STRINGS
      VARNAME = TEXT_LINE(1:EQLOC-1)
      VARVAL  = ADJUSTL(TEXT_LINE(EQLOC+1:LENGTH))
      LVARVAL = LEN_TRIM(VARVAL)
      IF(LVARVAL .eq. 0) CALL PERROR(6,'IN DATA PARAMETER FILE', 'VARIABLE LINE'//ERRSTRING//' HAS NO ASSOCIATED VALUE')

      ! DETERMINE TYPE OF VARVAL
      ! CHECK FOR LOGICAL
      IF(((VARVAL(1:1) .eq. "T") .OR. (VARVAL(1:1) .eq. "F")) .and. (LVARVAL .eq. 1))THEN 
        VARTYPE = "logical"
        IF(VARVAL(1:1) .eq. "T") LOGVAL = .TRUE.
        RETURN
      END IF

      ! CHECK IF IT IS A STRING  (CONTAINS CHARACTERS OTHER THAN 0-9,+,-,e,E,.)
      DO ii = 1, LVARVAL
        IF(INDEX(NUMCHARS,VARVAL(ii:ii)) == 0) VARTYPE = "string" 
      END DO


      ! PROCESS STRING (MAY BE MULTIPLE)
      IF(VARTYPE .eq. "string") THEN
        TSTRING = VARVAL
        STRINGVAL(1) = TSTRING 
        NVAL = 1
        ONFRAG = .TRUE.
        DO ii = 1, LVARVAL
          IF(VARVAL(ii:ii) .ne. " ")THEN
            FRAG(NVAL) = TRIM(FRAG(NVAL))//VARVAL(ii:ii)
            ONFRAG = .TRUE.
          ELSE
            IF(ONFRAG) NVAL = NVAL + 1
            ONFRAG = .FALSE.
          END IF
        END DO
        DO ii = 1, NVAL
          STRINGVAL(ii + 1) = TRIM(FRAG(ii))
        END DO
        RETURN
      END IF

      ! CHECK IF IT IS A FLOAT
      DOTLOC = INDEX(VARVAL,".")
      IF(DOTLOC .ne. 0) THEN
        VARTYPE = "float"
      ELSE
        VARTYPE = "integer"
      END IF
   
    ! FRAGMENT INTO STRINGS FOR MULTIPLE VALUES
      NP = 1
      ONFRAG = .TRUE.
      DO ii = 1, LVARVAL
        IF (VARVAL(ii:ii) .ne. " ")THEN 
          FRAG(NP) = TRIM(FRAG(NP))//VARVAL(ii:ii)
          ONFRAG = .TRUE.
        ELSE
          IF(ONFRAG) NP = NP + 1
          ONFRAG = .FALSE.
        END IF
      END DO

      ! EXTRACT NUMBER(S) FROM CHARACTER STRINGS
      NVAL = NP
      DO ii = 1, NP
        TEMP = TRIM(FRAG(ii))
        IF (VARTYPE .eq. "float") THEN 
          READ(TEMP, *) REALVAL(ii)
        ELSE
          READ(TEMP, *) INTVAL(ii)
        END IF
      END DO
  
    END SUBROUTINE GET_VAL 


  !==============================================================================|
  integer function SCAN_FILE(FNAME,VNAME,ISCAL,FSCAL,IVEC,FVEC,CVEC,NSZE,CVAL,LVAL)
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

  !   OPTIONAL (MUST PROVIDE ONE)
  !        ISCAL = INTEGER SCALAR
  !        FSCAL = FLOAT SCALAR
  !        CVAL = CHARACTER VARIABLE
  !        LVAL = LOGICAL VARIABLE
  !        IVEC = INTEGER VECTOR **
  !        FVEC = FLOAT VECTOR **
  !        CVEC = STRING VECTOR **
  !      **NSZE = ARRAY SIZE (MUST BE PROVIDED WITH IVEC/FVEC)

     USE MOD_PREC
     IMPLICIT NONE
     CHARACTER(LEN = *) :: FNAME, VNAME
     INTEGER, INTENT(INOUT), OPTIONAL :: ISCAL, IVEC(*)
     REAL(SP), INTENT(INOUT), OPTIONAL :: FSCAL, FVEC(*)
     CHARACTER(LEN=80), OPTIONAL :: CVAL, CVEC(*)
     LOGICAL, INTENT(INOUT), OPTIONAL :: LVAL
     INTEGER, INTENT(INOUT), OPTIONAL :: NSZE
   
     REAL(SP) REALVAL(150)
     INTEGER  INTVAL(150)
     CHARACTER(LEN=20 ) :: VARNAME
     CHARACTER(LEN=80 ) :: STRINGVAL(150),TITLE
     CHARACTER(LEN=80 ) :: INPLINE
     CHARACTER(LEN=400) :: TLINE
     CHARACTER(LEN=7  ) :: VARTYPE
     CHARACTER(LEN=20 ), DIMENSION(200) :: SET
     INTEGER :: I, NVAL, J, NSET, NLINE, NREP
     LOGICAL :: SETYES, ALLSET, CHECK, LOGVAL


     SCAN_FILE = 0

    ! OPEN THE INPUT FILE
    INQUIRE(FILE = TRIM(FNAME), EXIST = CHECK)
    IF (.NOT. CHECK) THEN
      SCAN_FILE = -1
      RETURN
    END IF

    OPEN(10, FILE=TRIM(FNAME)); REWIND(10) 


  ! SCAN THE FILE FOR THE VARIABLE NAME
     NSET = 0
     NLINE = 0
     DO WHILE (.TRUE.)
       TLINE(1:LEN(TLINE)) = ' ' 
       NREP  = 0
       NLINE = NLINE + 1
       READ(10,'(a)',END=20) INPLINE
       TLINE(1:80) = INPLINE(1:80)

  ! PROCESS LINE CONTINUATIONS
   110 CONTINUE
       I = LEN_TRIM(INPLINE)
       IF (I .ne. 0) THEN
  !#    if defined (COMPAQ) || defined (INTEL) || defined (IRIX)
       !IF( INPLINE(I-1:I) == '\\')THEN
  !#    else
       IF ( INPLINE(I-1:I) == '\\\\') THEN
  !#    endif
         NREP = NREP + 1
         READ(10,'(a)',END=20) INPLINE
         NLINE = NLINE + 1
         TLINE( NREP*80 + 1 : NREP*80 +80) = INPLINE(1:80)
         GOTO 110
       END IF
       END IF
       IF(NREP > 4)CALL PERROR(6,"CANNOT HAVE > 4 LINE CONTINUATIONS")

    ! REMOVE LINE CONTINUATION CHARACTER \\
       IF (NREP .gt. 0)THEN
         DO I = 2, LEN_TRIM(TLINE)
  !#    if defined (COMPAQ) || defined (INTEL) || defined (IRIX)
           !IF( TLINE(I-1:I) == '\\') TLINE(I-1:I) = '  '
  !#        else
           IF( TLINE(I-1:I) == '\\\\') TLINE(I-1:I) = '  '
  !#        endif
         END DO
       END IF
       
    ! PROCESS THE LINE
    CALL GET_VAL(NLINE, LEN_TRIM(TLINE), ADJUSTL(TLINE), VARNAME, VARTYPE, LOGVAL, STRINGVAL, REALVAL, INTVAL, NVAL)

  ! IF VARNAME MATCHES, PROCESS VARIABLE AND ERROR-CHECK

       IF(TRIM(VARNAME) .eq. TRIM(VNAME))THEN

         IF(PRESENT(ISCAL))THEN
           IF(VARTYPE .eq. 'integer')THEN
             ISCAL = INTVAL(1)
             RETURN
           ELSE
             SCAN_FILE = -3
           END IF
         ELSE IF(PRESENT(FSCAL))THEN
           IF(VARTYPE .eq. 'float')THEN
             FSCAL = REALVAL(1)
             RETURN
           ELSE
             SCAN_FILE = -3
           END IF
         ELSE IF(PRESENT(CVAL))THEN
           IF(VARTYPE .eq. 'string')THEN
             CVAL = STRINGVAL(1) 
             RETURN
           ELSE
             SCAN_FILE = -3
           END IF
         ELSE IF(PRESENT(LVAL))THEN
           IF(VARTYPE .eq. 'logical')THEN
             LVAL = LOGVAL 
             RETURN
           ELSE
             SCAN_FILE = -3
           END IF
         ELSE IF (PRESENT(IVEC)) THEN
           IF (NVAL .gt. 1) THEN
             IF (VARTYPE .eq. 'integer') THEN
               IVEC(1:NVAL) = INTVAL(1:NVAL) 
               NSZE = NVAL
               RETURN
             ELSE
               SCAN_FILE = -3
             END IF
             ELSE
             SCAN_FILE = -4 
           END IF
         ELSE IF (PRESENT(FVEC)) THEN
           IF (NVAL .gt. 1) THEN
             IF (VARTYPE .eq. 'float') THEN
               FVEC(1:NVAL) = REALVAL(1:NVAL) 
               NSZE = NVAL           
               RETURN
             ELSE
               SCAN_FILE = -3
             END IF
           ELSE
             SCAN_FILE = -4 
           END IF
         ELSE IF (PRESENT(CVEC)) THEN
           IF(NVAL .gt. 0)THEN
             IF (VARTYPE .eq. 'string') THEN
               CVEC(1:NVAL) = STRINGVAL(2:NVAL+1)
               NSZE = NVAL 
               RETURN
             ELSE
               SCAN_FILE = -3
             END IF
           ELSE
             SCAN_FILE = -4
           END IF
         ELSE
           SCAN_FILE = -5
         END IF
       END IF  !!VARIABLE IS CORRECT
            
     END DO !!LOOP OVER INPUT FILE
   20 CLOSE(10) 
     SCAN_FILE = -2

     RETURN 
     END FUNCTION SCAN_FILE


    SUBROUTINE PERROR(IOUT, ER1, ER2, ER3, ER4) !--Huang change -- add from FVCOM "mod_utils.F"

       IMPLICIT NONE
       INTEGER :: IOUT
       CHARACTER(LEN = *)             :: ER1
       CHARACTER(LEN = *), OPTIONAL   :: ER2
       CHARACTER(LEN = *), OPTIONAL   :: ER3
       CHARACTER(LEN = *), OPTIONAL   :: ER4

       WRITE(IOUT, *) '==================ERROR=================================='
       WRITE(IOUT, *) ER1
       IF (PRESENT(ER2)) WRITE(IOUT, *) ER2
       IF (PRESENT(ER3)) WRITE(IOUT, *) ER3
       IF (PRESENT(ER4)) WRITE(IOUT, *) ER4
       WRITE(IOUT, *) '========================================================='
       STOP
       RETURN

    END SUBROUTINE PERROR

END MODULE MOD_INP

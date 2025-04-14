module io

    use variables, only : sp
    implicit none

contains
    
    function scanInteger(file, key) result(value)
        character(len = *), intent(in) :: key, file
        character(len = 80) :: buffer
        integer :: value
        buffer = scanLines(file, key)
        read(buffer, *) value
    end function

    function scanReal(file, key) result(value)
        character(len = *), intent(in) :: key, file
        character(len = 80) :: buffer
        real(sp) :: value
        buffer = scanLines(file, key)
        read(buffer, *) value
    end function

    function scanString(file, key) result(value)
        character(len = *), intent(in) :: key, file
        character(len = 80) :: value
        value = scanLines(file, key)
    end function

    function scanLogical(file, key) result(value)
        character(len = *), intent(in) :: key, file
        character(len = 80) :: buffer
        logical :: value
        buffer = scanLines(file, key)
        value = merge(.true., .false., buffer(1:1) == "T")
    end function

    function scanLines(file, key) result(value)
        character(len = *), intent(in) :: key, file
        character(len = 80) :: value, line
        integer :: fid = 10, split, length, lines=0
        open(fid, file=trim(file))
        rewind(fid)
        do while (.true.)
            lines = lines + 1
            read(fid,'(a)', end=20) line
            length = len_trim(line)
            split = index(line, "=")
            if (key == trim(line(1:split-1))) then
                value = adjustl(line(split+1:length))
                close(fid)
                return
            end if
        end do
        20 close(fid)
    end function
end module

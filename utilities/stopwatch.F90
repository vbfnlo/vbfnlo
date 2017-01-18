module stopwatch
    integer starttime, stoptime, count_rate
    contains

subroutine starttimer()
  implicit none
  call system_clock(COUNT=starttime,COUNT_RATE=count_rate)
end subroutine

subroutine stoptimer()
  implicit none
  call system_clock(COUNT=stoptime)
end subroutine

subroutine printduration()
  implicit none
  integer totalsec, seconds, minutes, hours
  totalsec = (stoptime-starttime)/count_rate
  minutes = totalsec/60
  seconds = mod(totalsec, 60)
  hours = minutes/60
  minutes = mod(minutes, 60)
  print '("  Runtime: ",i0,":",i0.2,":",i0.2," (",i0," seconds)")', hours, minutes, seconds, totalsec
end subroutine

subroutine printtime(prefix)
  implicit none
  character(*), intent(in) :: prefix
  character(5)  :: zone
  integer,dimension(8) :: values
  ! using keyword arguments
  call date_and_time(ZONE=zone,VALUES=values)
  print '(a11,i4,"-",i0.2,"-",i0.2,"T",i0.2,":",i0.2,":",i0.2,a5)', &
       prefix, values(1:3),            values(5:7),           zone
end subroutine

endmodule

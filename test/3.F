      integer, parameter :: n = 3
      double precision x(n), y(n)
      x = (/ 1.0D0, 2.0D0, 3.0D0 /)
      y = (/ 5.0D0, 6.0D0, 7.0D0 /)
      open(3, file='fort.3', status='new', form='unformatted')
      write(3) x
      write(3) y
      close(3)
      end


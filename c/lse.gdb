file /path/to/liblse.so

# Define a function to create LSE
define create_lse
  set $ngauss = $arg0
  set $nmax = $arg1
  set $epsilon = $arg2
  set $lse = (LSE*)lse_malloc($ngauss, $nmax, $epsilon)
  print $lse
end

# Define other helper functions
define compute
  call lse_compute($lse, $arg0)
end

define run_tmat
  call lse_tmat($lse)
end

define cleanup
  call lse_free($lse)
  set $lse = 0
end

program SimulatedEx

  !compilaion: $ make -f makesiman

  !Running for a):

  !$ ./SimulatedEx.exe 
  !Seed from /dev/urandom: 282062966

  !Scheme a, b or c?
  !a

  !Name of datafile that has the points?
  !cities
  !How many points?
  !20

  !Similarly for b) and c), but without the data file question.

  !Solving travelling salesman problem with Monte Carlo Simulated Annealing
  !method. Salesman sells vacuums and travels by helicopter. Travel speed is
  ! 500km/hour and selling a vacuum takes 20 minutes (no need to sell vacuum in
  !home city). Code has three schemes, input (a) if you want to read in the data
  !file included with the cities, (b) in order to create a random set of points
  !in a unit square and (c) in order to create points randomly on the unit
  !circle. Output is distance travelled and code running time. Also for case a),
  !total time spent travelling is printed out and the order in which the cities 
  !where visited (total time includes sales of vacuums).

  !T_o is chosen to be 10000, in 10000000 steps the T is lowered linearly
  !to 0.1*T_o, and after this T_o is set as 0.8T_o and another 10000000 steps
  !are simulated. This is done 100 times. Arrays named 'indx' and 'distances'
  !contain the information of the sequence of points and their distances from
  !each other. Some total distances of the routes are printed into data files
  !for plotting.

  !Modules used
  use mtdefs
  use mtmod
  implicit none

  double precision, allocatable :: distances(:,:), points(:,:)
  double precision :: T, T_0, start, end, dx, dy, dist, dist_0, v, o, p, pi
  integer :: n, iarg, i, j, k, N_mc, N_siman, u, w, a, b
  integer, allocatable :: indx(:)
  character (len=15) :: arg, scheme, c, d
  character (len=15), allocatable :: cities(:)

  ! Initialize RNG
  call sgrnd(getseed(info=1))

  call cpu_time(start)

  ! initial T, and number of steps
  T_0 = 10000
  N_mc = 100
  N_siman = 10000000

  ! Read in which scheme a, b or c are we doing
  print *,
  print *, "Scheme a, b or c?"
  read(*,*) scheme
  if ((scheme/='a').and.(scheme/='b').and.(scheme/='c')) then
     print *, "Give a, b or c!"
     stop
  end if

!---------- If we have scheme 'a' ---------- 
  if (scheme=='a') then
     
     ! Read in how many points and file name
     print *,
     print *, "Name of datafile that has the points?"
     read(*,*) arg
     print *, "How many points?"
     read(5,*) n

     ! Allocate arrays
     allocate(points(n,2), distances(n,n), cities(n), indx(n))
     
     ! Create the index array where we suffle the cites
     do i = 1,n
        indx(i) = i
     end do

     ! Read in the file that has the points
     open(unit=1, file=arg, status='unknown')
     do i=1,n
        read(1,*) points(i,1:2), cities(i)
     end do
     close(1)

     ! Calculate the distances of the points into the distances-array
     do i=1,n
        do j=1,n
           dx = points(j,1)-points(i,1)
           dy = points(j,2)-points(i,2)
           distances(j,i) = sqrt(dx*dx+dy*dy)
        end do
     end do

     ! Calculating the initial distance
     dist_0 = 0
     do i = 1,n-1
        dist_0 = dist_0+distances(i,i+1)
     end do
     dist_0 = dist_0+distances(n,1) ! Last step

     ! Do the simulated annealing starting with the cities in the order
     ! they appear in the list
     open(unit=2, file='caseA', status='unknown')
     do i = 1,N_mc
        T = T_0
        do j = 1,N_siman
           
           ! Choose the two cities to exchange
           u = igrnd(1,n)
           w = igrnd(1,n)
           if (u==w) then ! If the same city is selected
              w = igrnd(1,n)
           end if
           ! Exchanging the cities
           a = indx(u)
           b = indx(w)
           indx(u) = b
           indx(w) = a
           
           ! Calculate new distance
           dist = 0
           do k = 1,n-1
              dist = dist+distances(indx(k),indx(k+1))
           end do
           dist = dist+distances(indx(n),indx(1)) ! Last step

           ! Deciding if we accept the change
           if ((dist-dist_0)<=0) then
              c = cities(u) 
              d = cities(w) 
              cities(u) = d
              cities(w) = c
              dist_0 = dist
           else if ((dist-dist_0)>0) then
              v = grnd()
              if (v<exp((dist_0-dist)/T)) then
                 c = cities(u) 
                 d = cities(w) 
                 cities(u) = d
                 cities(w) = c
                 dist_0 = dist
              else
                 indx(w) = b
                 indx(u) = a
              end if
           end if
           
           ! Write some values to a file
           if ((mod(j,10000)==1).or.(j==1)) then
              write(2,*) dist_0
           end if
           ! Lower the temperature
           T = T-0.00000001*T_0
           
        end do
        
        ! Raise the temperature
        T_0 = T_0*(0.8d+0)
     end do
     close(2)

     ! print out the sequence of cities
     do i=1,n
        print *,
        print *, cities(i)
     end do
     print *, "Distance travelled:"
     print *, dist_0
     print *, "Time taken:"
     print *, (500/dist_0)+0.3*19
     print *,

  end if


!---------- If we have scheme 'b' ----------

  if (scheme=='b') then

     ! Read in how many points and file name
     print *,
     print *, "How many points?"
     read(5,*) n
     
     allocate(points(n,2), distances(n,n), indx(n))

     ! Create the index array where we suffle the points
     do i = 1,n
        indx(i) = i
     end do


     ! Create random points
     do i=1,n
        o = grnd()
        p = grnd()
        points(i,1) = o
        points(i,2) = p
     end do
     
     ! Calculate the distances of the points into the distances-array
     do i=1,n
        do j=1,n
           dx = points(j,1)-points(i,1)
           dy = points(j,2)-points(i,2)
           distances(j,i) = sqrt(dx*dx+dy*dy)
        end do
     end do
  
     ! Calculating the initial distance
     dist_0 = 0
     do i = 1,n-1
        dist_0 = dist_0+distances(i,i+1)
     end do
     dist_0 = dist_0+distances(n,1) ! Last step

     ! Do the simulated annealing 
     open(unit=3, file='caseB', status='unknown')
     do i = 1,N_mc
        T = T_0
        do j = 1,N_siman
           
           ! Choose the two points to exchange
           u = igrnd(1,n)
           w = igrnd(1,n)
           if (u==w) then ! If the same point is selected
              w = igrnd(1,n)
           end if    

           ! Exchanging the points
           a = indx(u)
           b = indx(w)
           indx(u) = b
           indx(w) = a
           
           ! Calculate new distance
           dist = 0
           do k = 1,n-1
              dist = dist+distances(indx(k),indx(k+1))
           end do
           dist = dist+distances(indx(n),indx(1)) ! Last step
           
           ! Deciding if we accept the change
           if ((dist-dist_0)<=0) then
              dist_0 = dist
           else if ((dist-dist_0)>0) then
              v = grnd()
              if (v<exp((dist_0-dist)/T)) then
                 dist_0 = dist
              else
                 indx(w) = b
                 indx(u) = a
              end if
           end if
           
           ! Write some values to a file
           if ((mod(j,10000)==1).or.(j==1)) then
              write(3,*) dist_0
           end if
           ! Lower the temperature
           T = T-0.00000001*T_0
           
        end do

        ! Raise the temperature
        T_0 = T_0*(0.8d+0)
     end do
     close(3)

     ! print out distance 
     print *, "Distance travelled:"
     print *, dist_0
     print *,

  end if

!---------- If we have scheme 'c' ----------

  if (scheme=='c') then

     pi = 4.0*atan(1.0)
     
     ! Read in how many points and file name
     print *,
     print *, "How many points?"
     read(5,*) n
     
     allocate(points(n,2), distances(n,n), indx(n))
     
     ! Create the index array where we suffle the points
     do i = 1,n
        indx(i) = i
     end do
     
     
     ! Create random points on unit circle
     do i=1,n
        o = grnd()*2*pi
        points(i,1) = cos(o)
        points(i,2) = sin(o)
     end do
     
     ! Calculate the distances of the points into the distances-array
     do i=1,n
        do j=1,n
           dx = points(j,1)-points(i,1)
           dy = points(j,2)-points(i,2)
           distances(j,i) = sqrt(dx*dx+dy*dy)
        end do
     end do
     
     ! Calculating the initial distance
     dist_0 = 0
     do i = 1,n-1
        dist_0 = dist_0+distances(i,i+1)
     end do
     dist_0 = dist_0+distances(n,1) ! Last step

     ! Do the simulated annealing 
     open(unit=4, file='caseC', status='unknown')
     do i = 1,N_mc
        T = T_0
        do j = 1,N_siman
           
           ! Choose the two points to exchange
           u = igrnd(1,n)
           w = igrnd(1,n)
           if (u==w) then ! If the same point is selected
              w = igrnd(1,n)
           end if    

           ! Exchanging the points
           a = indx(u)
           b = indx(w)
           indx(u) = b
           indx(w) = a
           
           ! Calculate new distance
           dist = 0
           do k = 1,n-1
              dist = dist+distances(indx(k),indx(k+1))
           end do
           dist = dist+distances(indx(n),indx(1)) ! Last step
           
           ! Deciding if we accept the change
           if ((dist-dist_0)<=0) then
              dist_0 = dist
           else if ((dist-dist_0)>0) then
              v = grnd()
              if (v<exp((dist_0-dist)/T)) then
                 dist_0 = dist
              else
                 indx(w) = b
                 indx(u) = a
              end if
           end if

          ! Write some values to a file
           if ((mod(j,10000)==1).or.(j==1)) then
              write(4,*) dist_0
           end if
           ! Lower the temperature
           T = T-0.00000001*T_0
           
        end do

        ! Raise the temperature
        T_0 = T_0*(0.8d+0)
     end do
     close(4)

     ! print out distance 
     print *, "Distance travelled:"
     print *, dist_0
     print *,

  end if 


  call cpu_time(end)
  print *, (end-start)/60

  


end program SimulatedEx
  
  
  

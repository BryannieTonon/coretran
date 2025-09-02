program geographic_network_kdtree
  use variableKind, only: i32, r64
  use m_allocate, only: allocate
  use m_deallocate, only: deallocate
  use m_random, only: rngUniform
  use m_KdTree, only: KdTree, KdTreeSearch
  use dArgDynamicArray_class, only: dArgDynamicArray
  use rdn_mod
  !use m_KdTree, only: KdTree, KdTreeSearch
  !use dArgDynamicArray_Class, only: dArgDynamicArray

  implicit none

  ! COMANDO PARA COMPILAR - gfortran -I./include app/mod_gerador.f90 app/kdTree_simple.f90 ./lib/libcoretran.so -o kdTree_simple -Wl,-rpath=./lib
  ! ./kdTree_simple

  integer(i32), parameter :: N = 30 ! Num de vértices
  real(r64), parameter :: L = 10.0_r64 ! Tamanho do quadrado
  real(r64), parameter :: rmin = 0.1_real64
  real(r64), parameter :: rmax = L / 5.0_real64
  real(r64) ::  R ! Raio
  integer(i32) :: i, j, k, unidade_arquivo
  character(len=*), parameter :: formatador = '(*(g0,x))'
  character(len=*), parameter :: formatador2 = '(*(g0,1X))'


  real(r64), allocatable :: x(:), y(:), raio(:)
  integer(i32), allocatable :: degree(:)
  integer(i32), allocatable :: neighborList(:,:)

  type(KdTree) :: tree
  type(dArgDynamicArray) :: neighbors, da
  real(r64) :: dx, dy, dist2

  type(KdTreeSearch) :: search

  ! Alocação e sorteio
  
  call allocate(x, N) ! função do coretran
  call allocate(y, N)
  allocate(degree(N), raio(N))
  allocate(neighborList(N,N))
  neighborList = 0
  degree = 0

  ! Sorteios - PERGUNTA - como trocar todo esse processo de gerar os rdn?
  call generate_random_geo(N, L, rmin, rmax, x, y, raio)

  ! Constrói KD-Tree
  tree = KdTree(x, y)

!!! K-Nearest to (0, 0), for 3D add z, and zQuery
!!da = search%kNearest(tree, x, y, [z], xQuery = 0.d0, yQuery = 0.d0, [zQuery = 0.d0], k = 10) 
!!! Search for all points within a given distance
!!da = search%kNearest(tree, x, y, [z], xQuery = 0.d0, yQuery = 0.d0, [zQuery = 0.d0], radius = 10.d0)
    open(unit=unidade_arquivo,file='listas-de-pares-coretran.dat',action='write',status='unknown')

  ! Busca vizinhos dentro dentro de um raio R
  do i = 1, N ! Pega todos os pontos próximos (k=N-1) e filtra pelo raio
     R = raio(i)
     neighbors = search%kNearest(tree, x, y, xQuery = x(i), yQuery = y(i), radius = R)
     k = 0
     do j = 1, neighbors%size()
        dx = x(neighbors%i%values(j)) - x(i)
        dy = y(neighbors%i%values(j)) - y(i)
        dist2 = dx*dx + dy*dy
        ! Esse if aborda auto-conexão, mas não multiplas conexões - PERGUNTA - devido aos diferentes raios, como fica as multiplas conexoes?
        if ( (dist2 <= R*R) .and. (i /= neighbors%i%values(j))) then
            k = k + 1
            neighborList(i,k) = neighbors%i%values(j)
            write(unidade_arquivo,formatador) i, neighbors%i%values(j)
        end if
     end do
     degree(i) = k
    end do
    close(unidade_arquivo)

    open(unit=unidade_arquivo,file='coordenadas-raio-coretran.dat',action='write',status='unknown')
    do i=1 , N
        write(unidade_arquivo,formatador2) degree(i), x(i), y(i), raio(i)
    end do
    close(unidade_arquivo)

  call deallocate(x)
  call deallocate(y)
  deallocate(degree)
  deallocate(neighborList)
  call tree%deallocate()

end program geographic_network_kdtree

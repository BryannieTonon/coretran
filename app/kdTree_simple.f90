program geographic_network_kdtree
  use variableKind, only: i32, r64
  use m_allocate, only: allocate
  use m_deallocate, only: deallocate
  use m_random, only: rngUniform
  use m_KdTree, only: KdTree, KdTreeSearch
  use dArgDynamicArray_class, only: dArgDynamicArray
  use rndgen_mod
  use read_tools
 
  implicit none

  ! COMANDO PARA COMPILAR - gfortran -I./include app/rndgen.f90 app/read_tools.f90 app/kdTree_simple.f90 ./lib/libcoretran.so -o kdTree_simple -Wl,-rpath=./lib
  ! ./kdTree_simple

  !integer(i32), parameter :: N = 1000 ! Num de vértices
  !real(r64), parameter :: L = 500.0_r64 ! Tamanho do quadrado
  real(r64), parameter :: fator_escala = 10.0_r64

  real(r64), parameter :: rmin = 0.1_r64
  !real(r64), parameter :: rmax = L / fator_escala
  real(r64) ::  R ! Raio
  character(len=*), parameter :: formatador = '(*(g0,x))'
  integer(i32) :: node_1, node_2, N
  real(r64) :: dx, dy, dist2, L, rmax
  integer(i32) :: i, j, k, unidade_arquivo, int_fatorescala, int_L
  real(r64), allocatable :: x(:), y(:), raio(:)
  integer(i32), allocatable :: degree(:)
  integer(i32), allocatable :: neighborList(:,:)
character(len=200) :: comando
  integer :: seed = 294727492  
  type(KdTree) :: tree
  type(dArgDynamicArray) :: neighbors, da
  type(rndgen) :: rnd
  type(KdTreeSearch) :: search

  ! INPUT
    if (command_argument_count() /= 2) then
        print*, 'Faltam argumentos:'
        print*,"for N in {100..500..10}; do L=15"
        print*, "./kdTree_simple $N $L done"
        stop
    end if

    L = read_int(1)
    N = read_int(2)
    rmax = L / fator_escala

! Alocação e sorteio
  call rnd%init(seed)
  call allocate(x, N) 
  call allocate(y, N)
  allocate(degree(N), raio(N))
  allocate(neighborList(N,N))
  neighborList = 0
  degree = 0
  int_fatorescala = int(fator_escala)
  int_L = int(L)
! Sorteios 
!generator%rnd_array(n) for a real array with numbers in range [0,1) with size n
!generator%rnd_array(n,i1,i2) for a integer array with numbers in the range [i1, i2] with size n
!generator%rnd_array(n,r1,r2) for a real array with numbers in the range [r1, r2) with size n

! RAIO
  raio = rnd%rnd_array(N,rmin,rmax)
! COORDENADAS
x = rnd%rnd_array(N,0.0_r64,L)
y = rnd%rnd_array(N,0.0_r64,L)

! Constrói KD-Tree
tree = KdTree(x, y)

!!! K-Nearest to (0, 0), for 3D add z, and zQuery
!!da = search%kNearest(tree, x, y, [z], xQuery = 0.d0, yQuery = 0.d0, [zQuery = 0.d0], k = 10) 
!!! Search for all points within a given distance
!!da = search%kNearest(tree, x, y, [z], xQuery = 0.d0, yQuery = 0.d0, [zQuery = 0.d0], radius = 10.d0)

open(newunit=unidade_arquivo,file='listas-de-pares-coretran.dat',action='write',status='unknown')

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
            node_1 = i
            node_2 = neighbors%i%values(j)
            call swap(node_1, node_2)
            write(unidade_arquivo,formatador) node_1, node_2
        end if
    end do
    degree(i) = k
    end do
    close(unidade_arquivo)

    ! Coordendas e raios para plotar a rede
    !open(newunit=unidade_arquivo,file='coordenadas-raio-coretran.dat',action='write',status='unknown')
    !do i=1 , N
    !    write(unidade_arquivo,formatador) degree(i), x(i), y(i), raio(i)
    !end do
    !close(unidade_arquivo)


    !comando = 'sort -n -k1,1 -k2,2 Nfixo_500-L_100-2-uniq.dat | uniq > Nfixo_500-L_100-2-TESTE.dat'
    comando = 'sort -n -k1,1 -k2,2 listas-de-pares-coretran.dat | uniq > ' // &
          'Lfixo_'// aux_name_int(int_L)// &
          '-N_'// aux_name_int(N)// &
          '-' //aux_name_int(int_fatorescala) // '.dat'

    call system(trim(comando))

    ! esse aqui deu errado
    !call system('sort -n '//'listas-de-pares-coretran.dat'//' | uniq > '//'Lfixo_'//aux_name_int(L)//'-N_'//aux_name_int(N)//'-'//aux_name_int(int_fatorescala)//'-uniq.dat')
    !call system( &
   !'sort -n listas-de-pares-coretran.dat | uniq > ' // &
   !'Nfixo_'// aux_name_int(N) // &
   !'-L_'// aux_name_int(int_L) // &
   !'-'// aux_name_int(int_fatorescala) // &
   !'-arquivo.dat' )

  call deallocate(x)
  call deallocate(y)
  deallocate(degree)
  deallocate(neighborList)
  call tree%deallocate()

contains

    subroutine swap(a, b)
        integer(i32), intent(inout) :: a, b
        integer(i32) :: temp

        if (a > b) then
            temp = a
            a = b
            b = temp
        end if
    end subroutine swap

end program geographic_network_kdtree

program geographic_network_kdtree
use variableKind, only: i32, r64
use m_allocate, only: allocate
use m_deallocate, only: deallocate
use m_random, only: rngUniform
use m_KdTree, only: KdTree, KdTreeSearch
use dArgDynamicArray_class, only: dArgDynamicArray
use rndgen_mod
use read_tools
!use tictoc_mod

implicit none

  ! COMANDOs PARA COMPILAR 

    ! gfortran -I./include app/rndgen.f90 app/read_tools.f90 app/kdTree_simple.f90 ./lib/libcoretran.so -o kdTree_simple -Wl,-rpath=./lib
    ! gfortran -I./include app/rndgen.f90 app/tictoc_mod.f90  app/read_tools.f90 app/kdTree_simple.f90 ./lib/libcoretran.so -o kdTree_simple -Wl,-rpath=./lib  
    ! ./kdTree_simple

    ! Mantendo um fixo e o outro variando
        ! for N in {100..700..20}; do     L=500;     ./kdTree_simple $L $N; done

    ! Variando N e L
        ! for N in {100..700..20}; do for L in {100..1000..100}; do ./kdTree_simple $L $N done done

    ! Variando N, L e fator_escala
        ! for N in {100..10000..20}; do for L in {100..1000..100}; do ./kdTree_simple $L $N done done
        ! for N in {5500..10000..500}; do     for L in {500..1000..500}; do         for F in {5..20..5}; do             ./kdTree_simple $L $N $F;         done;     done; done

  
  ! Parâmetros
  !integer(i32), parameter :: N = 1000 ! Num de vértices
  !real(r64), parameter :: L = 500.0_r64 ! Tamanho do quadrado
  !real(r64), parameter :: fator_escala = 2.0_r64

  real(r64), parameter :: rmin = 0.1_r64
  !real(r64), parameter :: rmax = L / fator_escala
  real(r64) ::  R ! Raio
  character(len=*), parameter :: formatador = '(*(g0,x))'
  integer(i32) :: node_1, node_2, N, fator_escala, n_arestas
  real(r64) :: dx, dy, dist2, L, rmax
  integer(i32) :: i, j, k, unidade_arquivo, unidade_tempo, int_L, int_fatorescala, num, id_amostra
  real(r64), allocatable :: x(:), y(:), raio(:)
  integer(i32), allocatable :: degree(:)
  integer(i32), allocatable :: neighborList(:,:) !, arestas(:,:), arestas_temp(:,:)
  character(len=200) :: comando
  integer :: seed 

  type(KdTree) :: tree
  type(dArgDynamicArray) :: neighbors, da
  type(rndgen) :: rnd
  type(KdTreeSearch) :: search
  
  !type(tictoc_t)    :: ctimer
  
  ! INPUT
if (command_argument_count() /= 4) then
    print*, 'Faltam argumentos:'
    print*,"for N in {100..500..10}; do L=15"
    print*, "./kdTree_simple $N $L done"
    stop
end if

L = read_int(1)
N = read_int(2)
fator_escala = read_int(3)
id_amostra = read_int(4)
!fator_escala = 4

rmax = L / (1.0_dp*fator_escala)

seed = 294727492 + id_amostra
    
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

! RAIO
  raio = rnd%rnd_array(N,rmin,rmax)
  
! COORDENADAS
    x = rnd%rnd_array(N,0.0_r64,L)
    y = rnd%rnd_array(N,0.0_r64,L)

! Constrói KD-Tree
    tree = KdTree(x, y)


open(newunit=unidade_arquivo,file='lista-de-pares.dat',action='write',status='unknown')
!open(newunit=unidade_tempo, file="tempo-CPU-L_"// aux_name_int(int_L)// &
!         "-N_"// aux_name_int(N)// "-"//aux_name_int(int_fatorescala) // &
!         "kdtree.dat",position="append",action='write',status='unknown')


!call ctimer%reset()
!call ctimer%tic()

do i = 1, N ! Pega todos os pontos próximos (k=N-1) e filtra pelo raio
    R = raio(i)
    ! vizinhos ordenados pela distância
    neighbors = search%kNearest(tree, x, y, xQuery = x(i), yQuery = y(i), radius = R)
    k = 0
    do j = 1, neighbors%size()
        dx = x(neighbors%i%values(j)) - x(i)
        dy = y(neighbors%i%values(j)) - y(i)
        dist2 = dx*dx + dy*dy

        ! REGRAS DE CONECTIVIDADE
        ! 1) Raio Fixo

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
!call ctimer%toc()
close(unidade_arquivo)
!write(unidade_tempo, formatador) N, int(L), fator_escala, ctimer%t_tot
!close(unidade_tempo)


! SAÍDAS

! Coordendas e raios para plotar a rede
! open(newunit=unidade_arquivo,file='coordenadas-raio-simples-teste.dat',action='write',status='unknown')
! do i=1 , N
!     write(unidade_arquivo,formatador) degree(i), x(i), y(i), raio(i)
! end do
! close(unidade_arquivo)

! lista de aresta
!comando = 'sort -n -k1,1 -k2,2 lista-de-pares.dat | uniq > lista-de-pares-simples-metricas.dat'
 comando = 'sort -n -k1,1 -k2,2 lista-de-pares.dat| uniq > ' // &
         'L_'// aux_name_int(int_L)// &
         '-N_'// aux_name_int(N)// &
         '-' //aux_name_int(int_fatorescala) //'-' //aux_name_int(id_amostra) //'-simples.dat'
call system(trim(comando))


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

    subroutine ordena_vertices(e, n) 
        integer(i32), intent(inout) :: e(n,2)
        integer(i32), intent(in) :: n
        integer(i32) :: i,j
        integer(i32) :: t1,t2

        do i = 1, n-1
            do j = i+1, n
                ! se o elemento da linha i coluna 1 for maior que o elemento da linha j coluna 1 ou
                ! se eles forem iguais e  elemento da linha i da coluna 2 for maior que o elemento da linha j da coluna 2
                if (e(i,1) > e(j,1) .or. (e(i,1) == e(j,1) .and. e(i,2) > e(j,2))) then
                    ! atribui o valor das variávies auxiliares
                    t1 = e(i,1)
                    t2 = e(i,2)
                    ! como e(i,1) é maior que e(j,1), inverte eles, matendo a correspondêcia com a segunda coluna
                    e(i,1) = e(j,1)
                    e(i,2) = e(j,2)
                    e(j,1) = t1
                    e(j,2) = t2
                end if
            end do
        end do
    end subroutine ordena_vertices

    subroutine escrever_lista(e, n, arquivo)
        integer(i4), intent(in) :: e(n,2), n, arquivo
        integer(i4) :: i
        character(len=*), parameter :: formatador = '(*(g0,x))'

        write(arquivo,formatador) e(1,1), e(1,2)
        do i = 2, n
            if (e(i,1) /= e(i-1,1) .or. e(i,2) /= e(i-1,2)) then
                write(arquivo,formatador) e(i,1), e(i,2)
            end if
        end do
    end subroutine escrever_lista

end program geographic_network_kdtree

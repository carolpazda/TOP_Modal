using BMesh, LMesh, LinearSolve, LFEM, LinearAlgebra

include("coluna2D.jl")

function Roda_coluninha()

    # Cria a malha
    nx = 2
    ny = 4
    mesh = Coluna2D(nx,ny,force=800E3)

    # Define os expoentes
    penal = 2.9
    q = 2.6

    ne = Get_ne(mesh)
    @show ne

    # Cria as variáveis de projeto e coloca "1" nas de interesse
    x = 0.5*ones(ne)
    x[1] = 1.0
    x[2] = 1.0
    x[9] = 1.0
    x[10] = 1.0
    x[11] = 1.0
    x[13] = 1.0
    x[14] = 1.0
    x[16] = 1.0
    x[17] = 1.0
    x[19] = 1.0
    x[20] = 1.0
    x[22] = 1.0

    #@show x

    # Parametrização da rigidez e da tensão
    kparam(xe::Float64,p=penal)=xe^p
    sparam(xe::Float64,p=penal,q=q)=xe^(p-q)

    # Solve the linear static equilibrium
    U, F, linsolve = Solve_linear(mesh,x,kparam)

    # Array with stresses - Default is at the center of the element
    sigma = Stresses(mesh,U,x,sparam)

    # Free dofs
    free_dofs = mesh.free_dofs[1]

    # Dimensão
    dim = Get_dim(mesh)
    ngls = dim*Get_nn(mesh)

    # Matrizes globais
    K = Global_K(mesh,x,kparam)
    Ks = Global_Ks(mesh,sigma)

    # Gls que interessam
    Kn = Array(K[free_dofs,free_dofs])
    Ksn = Array(Ks[free_dofs,free_dofs])

    modos = eigen(Kn,-Ksn)

    # Preparando para ver no gmsh
    ne = Get_ne(mesh)
    saida = "testando.pos"
    Gmsh_init(saida,mesh)
    Gmsh_element_scalar(mesh,x,saida,"Design variables")

    # Definindo o 'X', que no caso é o nosso PHI
    X = zeros(ngls,size(modos.vectors,2))

    for i=1:size(modos.vectors,2)
        vtemp = zeros(ngls)
        Expand_vector!(vtemp,modos.vectors[:,i],free_dofs)
        X[:,i] .= vtemp
        Gmsh_nodal_vector(mesh,X[:,i],saida,"Modo $(modos.values[i]) ")
    end

    return modos

end
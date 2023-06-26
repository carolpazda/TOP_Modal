using BMesh, LMesh, LinearSolve, LFEM, LinearAlgebra

function Tres_barras()

    etype = :truss2D
    nn = 4
    ne = 3
    coord = [0.0 0.0 ;
             1.0 0.0 ;
            2.0 0.0 ;
            1.0 3.0]

    connect = [1 4;
               2 4;
               3 4]

    bm = Bmesh2D(etype,nn,ne,coord,connect,1.0,1.0,1,1)

    nbc = [4 2 -1000.0]
    hebc = [1 1;
            1 2;
            2 1;
            2 2;
            3 1;
            3 2]

    mat = [Material(Ex=7E9,density=1.0,νxy=0.3,limit_stress=1E9)]
    geom = [Geometry(A=1E-4, thickness=1.0)]

    options = Dict{Symbol,Matrix{Float64}}()

    Mesh2D(bm,mat,geom,hebc,nbc,options=options)

end

function Roda_tres_barras()

    # Cria a malha
    mesh = Tres_barras()

    # Define os expoentes
    penal = 2.9
    q = 2.6

    # Cria as variáveis de projeto
    x = [1.0, 0.5, 1.0]

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
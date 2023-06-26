function Coluna2D(nx::Int64=2,ny::Int64=6,etype=:truss2D;   
     Lx=2.0, Ly=4.0,  force=100E3, A=0.0 ,Ex=210E9,
     νxy=0.0, density=7850.0,thickness=0.1)


    @assert iseven(nx) "Coluna2D: nx deve ser um número par" 

    @assert etype==:truss2D || etype==:solid2D "Coluna2D:: $etype must be :truss2D or :solid2D"

    # Generate the mesh
    if etype==:truss2D
        bmesh = Bmesh_truss_2D(Lx,nx,Ly,ny)
    else
        bmesh = Bmesh_solid_2D(Lx,nx,Ly,ny)
    end

    # Apoios nas extremidades de baixo
    hebc = [1 1;
            1 2;
            nx+1 1; 
            nx+1 2 ]
              

    # Generate the load information
    f = force
    no_f = (nx+1)*(ny) + round(Int64,nx/2) + 1
    nbc = [no_f 2 -f]

    @show nx, ny, nx+1, no_f

    # Todas as barras com seção circular e 2.5 cm de raio
    r = 2.5E-2
    A = pi*r^2

    # Material e geometria
    mat = [Material(Ex=Ex)]
    geom = [Geometry(A=A)]

    # Gera a malha e devolve um tipo Mesh2D
    Mesh2D(bmesh,mat,geom,hebc,nbc)

end
#
# Problema de teste para o artigo do CILAMCE
#
# Para A = pi(2.5E-2)^2
#      E = 210E9
#      F_tot = 100E3
#
#
# U[3 * (5 - 1) + 3] = -0.0008429435907309232
# max_s / 1.0e6 = 7.612102074464771
# min_s / 1.0e6 = -19.024231700498344
# minimum(P) = 1.0794879813691488e8
#
#
# 

function Bracket3D(nx::Int64=3,ny::Int64=2,nz::Int64=6,etype=:truss3D;   
                  Lx=3.0, Ly=2.0, Lz=6.0, force=100E3, A=0.0 ,Ex=210E9,
                   νxy=0.0, density=7850.0,thickness=0.1)



    @assert etype==:truss3D || etype==:solid3D "Bracket3D:: $etype must be :truss3D or :solid3D"

    # Generate the mesh
    if etype==:truss3D
        bmesh = Bmesh_truss_3D(Lx,nx,Ly,ny,Lz,nz)
    else
        bmesh = Bmesh_solid_3D(Lx,nx,Ly,ny,Lz,nz)
    end

    # Generate the supports: one at each corner of
    # plane Z=0
    ebc = [ 1 1 0.0;
            1 2 0.0;
            1 3 0.0;
            9 1 0.0;
            9 2 0.0;
            9 3 0.0;
           73 1 0.0;
           73 2 0.0;
           73 3 0.0;
           81 1 0.0;
           81 2 0.0;
           81 3 0.0]
    
    # Generate the load information
    f = force
    nbc = [44 3 -f ]

    # Todas as barras com seção circular e 2.5 cm de raio
    r = 2.5E-2
    A = pi*r^2

    # Material e geometria
    mat = [Material(Ex=Ex)]
    geom = [Geometry(A=A)]

    # Gera a malha e devolve um tipo Mesh3D
    Mesh3D(bmesh,mat,geom,Int.(ebc),nbc)

end

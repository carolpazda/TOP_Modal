
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

function Bridge3D(nx::Int64=8,ny::Int64=1,nz::Int64=1,etype=:truss3D;   
                  Lx=8.0, Ly=1.0, Lz=1.0, force=100E3, A=0.0 ,Ex=210E9,
                   νxy=0.0, density=7850.0,thickness=0.1)



    @assert etype==:truss3D || etype==:solid3D "Bridge:: $etype must be :truss3D or :solid3D"

    # Generate the mesh
    if etype==:truss3D
        bmesh = Bmesh_truss_3D(Lx,nx,Ly,ny,Lz,nz)
    else
        bmesh = Bmesh_solid_3D(Lx,nx,Ly,ny,Lz,nz)
    end

    # Generate the supports: one at each corner of
    # plane Z=0
    ebc = [1 1 0.0;
            1 2 0.0;
            1 3 0.0;
            nx+1 1 0.0;
            nx+1 2 0.0;
            nx+1 3 0.0;
            (nx+1)*ny+1 1 0.0; 
            (nx+1)*ny+1 2 0.0; 
            (nx+1)*ny+1 3 0.0 ;
            (nx+1)*(ny+1) 1 0.0 ;
            (nx+1)*(ny+1) 2 0.0 ;
            (nx+1)*(ny+1) 3 0.0 ]
                            
    
    # Generate the load information
    f = force/14
    nbc = [2 3 -f ; 
           3 3 -f ;
           4 3 -f ;
           5 3 -f ; 
           6 3 -f ; 
           7 3 -f ; 
           8 3 -f ; 
           11 3 -f ;
           12 3 -f ;
           13 3 -f ;
           14 3 -f ;
           15 3 -f ;
           16 3 -f ;
           17 3 -f ;
           20 2 -f ; 
           26 2 -f]

    # Todas as barras com seção circular e 2.5 cm de raio
    r = 2.5E-2
    A = pi*r^2

    # Material e geometria
    mat = [Material(Ex=Ex)]
    geom = [Geometry(A=A)]

    # Gera a malha e devolve um tipo Mesh3D
    Mesh3D(bmesh,mat,geom,Int.(ebc),nbc)

end

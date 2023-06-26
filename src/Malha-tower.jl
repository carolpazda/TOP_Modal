function Tower3D(nx::Int64=2,ny::Int64=2,nz::Int64=6,etype=:truss3D;   
                  Lx=2.0, Ly=2.0, Lz=4.0, force=1.0, A=0.0 ,Ex=210E9,
                   νxy=0.0, density=7850.0,thickness=0.1)



    @assert etype==:truss3D || etype==:solid3D "Exemplo3D:: $etype must be :truss3D or :solid3D"

    # Generate the mesh
    if etype==:truss3D
        bmesh = Bmesh_truss_3D(Lx,nx,Ly,ny,Lz,nz)
    else
        bmesh = Bmesh_solid_3D(Lx,nx,Ly,ny,Lz,nz)
    end

    # Generate the supports: one at each corner of
    # plane Z=0
    ebc =  [1 1 0.0; # Frente - primeiro
            1 2 0.0; # Frente - primeiro
            1 3 0.0; # Frente - primeiro
            ((nx/2)+1) 1 0.0; # Frente - segundo
            ((nx/2)+1) 2 0.0; # Frente - segundo
            ((nx/2)+1) 3 0.0; # Frente - segundo
            nx+1 1 0.0; # Frente - terceiro
            nx+1 2 0.0; # Frente - terceiro
            nx+1 3 0.0; # Frente - terceiro
            (nx*ny) 1 0.0; # Meio - primeiro
            (nx*ny) 2 0.0; # Meio - primeiro
            (nx*ny) 3 0.0; # Meio - primeiro
            ((nx*ny)+1) 1 0.0; # Meio - segundo
            ((nx*ny)+1) 2 0.0; # Meio - segundo
            ((nx*ny)+1) 3 0.0; # Meio - segundo
            ((nx+1)*ny) 1 0.0; # Meio - terceiro
            ((nx+1)*ny) 2 0.0; # Meio - terceiro
            ((nx+1)*ny) 3 0.0; # Meio - terceiro
            (nx+1)*ny+1 1 0.0; # Fundo - primeiro
            (nx+1)*ny+1 2 0.0; # Fundo - primeiro
            (nx+1)*ny+1 3 0.0; # Fundo - primeiro
            (nx+1)*(ny+1)-1 1 0.0; # Fundo - segundo
            (nx+1)*(ny+1)-1 2 0.0; # Fundo - segundo
            (nx+1)*(ny+1)-1 3 0.0; # Fundo - segundo
            (nx+1)*(ny+1) 1 0.0; # Fundo - terceiro
            (nx+1)*(ny+1) 2 0.0; # Fundo - terceiro
            (nx+1)*(ny+1) 3 0.0] # Fundo - terceiro
             
    # Generate the load information
    nbc = [41 3 -force]

    # Todas as barras com seção circular e 2.5 cm de raio
    r = 2.5E-2
    A = pi*r^2

    # Material e geometria
    mat = [Material(Ex=Ex)]
    geom = [Geometry(A=A)]

    # Gera a malha e devolve um tipo Mesh3D
    Mesh3D(bmesh,mat,geom,Int.(ebc),nbc)

end

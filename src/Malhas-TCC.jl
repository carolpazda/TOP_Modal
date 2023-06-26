#
# MALHAS UTILIZADAS NO TCC
#
#
#--------------------------------------------------------------------------------------------------
#
# MALHA 1 - BRIDGE 
#  < _ ____________ _>
#     |\/|\/|\/|\/|
#      
function Malha1(nx::Int64=8,ny::Int64=2,nz::Int64=1,etype=:truss3D;   
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
    nbc = [37 1 -f;
           45 1 f]

    # Todas as barras com seção circular e 2.5 cm de raio
    r = 2.5E-2
    A = pi*r^2

    # Material e geometria
    mat = [Material(Ex=Ex)]
    geom = [Geometry(A=A)]

    # Gera a malha e devolve um tipo Mesh3D
    Mesh3D(bmesh,mat,geom,Int.(ebc),nbc)

end

#
#  MALHA 2  
#  
# >---------------
#   \   / | \   /|
#    \ /  |  \ / |
#    / \  |  / \ |
#   /   \ | /   \|
# >---------------
#               |
#               \/
#                P
#
function Malha2(nx::Int64=8,ny::Int64=8,etype=:truss2D;   
    Lx=4.0, Ly=4.0,  force=100E3, A=0.0 ,Ex=210E9,
     νxy=0.0, density=7850.0,thickness=0.1)



    @assert etype==:truss2D || etype==:solid2D "Bracket:: $etype must be :truss2D or :solid2D"

    # Generate the mesh
    if etype==:truss2D
        bmesh = Bmesh_truss_2D(Lx,nx,Ly,ny)
    else
        bmesh = Bmesh_solid_2D(Lx,nx,Ly,ny)
    end

    # Generate x=0, y=0 / Ly ((nx+1)*ny+1)
    ebc = [1 1;
           1 2;
           (nx+1)*ny+1  1; 
           (nx+1)*ny+1  2]
                

    # Generate the load information
    f = force
    nbc = [nx+1 2 -f]

    # Todas as barras com seção circular e 2.5 cm de raio
    r = 2.5E-2
    A = pi*r^2

    # Material e geometria
    mat = [Material(Ex=Ex)]
    geom = [Geometry(A=A)]

    # Gera a malha e devolve um tipo Mesh2D
    Mesh2D(bmesh,mat,geom,ebc,nbc)

end

#
#  MALHA 3
#
#        |
#       \/
#        p
#     ________
#    |\ /|\ /|
#    |/ \|/ \|
#    |-------|
#    |\ /|\ /|
#    |/ \|/ \|
#   -----------
#   ///////////
function Malha3(nx::Int64=8,ny::Int64=8,etype=:truss2D;   
    Lx=2.0, Ly=4.0,  force=100E3, A=0.0 ,Ex=210E9,
    νxy=0.0, density=7850.0, thickness=0.1)


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
           2 1; 
           2 2;
           3 1;
           3 2;
           4 1;
           4 2;
           5 1;
           5 2;
           6 1;
           6 2;
           7 1;
           7 2;
           8 1;
           8 2;
           9 1;
           9 2]
             

   # Generate the load information
   f = force
   nbc = [77 2 -f]

   # Todas as barras com seção circular e 2.5 cm de raio
   r = 2.5E-2
   A = pi*r^2

   # Material e geometria
   mat = [Material(Ex=Ex)]
   geom = [Geometry(A=A)]

   # Gera a malha e devolve um tipo Mesh2D
   Mesh2D(bmesh,mat,geom,hebc,nbc)

end

#
#  MALHA 4
#
#        
# p-> ________
#    |\ /|\ /|
#    |/ \|/ \|
#    |-------|
#    |\ /|\ /|
#    |/ \|/ \|
#   -----------
#   ///////////
function Malha4(nx::Int64=2,ny::Int64=2,etype=:truss2D;   
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
           2 1; 
           2 2;
           3 1;
           3 2]
             

   # Generate the load information
   f = force
   nbc = [7 1 f]

   # Todas as barras com seção circular e 2.5 cm de raio
   r = 2.5E-2
   A = pi*r^2

   # Material e geometria
   mat = [Material(Ex=Ex)]
   geom = [Geometry(A=A)]

   # Gera a malha e devolve um tipo Mesh2D
   Mesh2D(bmesh,mat,geom,hebc,nbc)

end


#
#  MALHA 5
#
#             |
#             \/
#             p
#      ________________  
#    /| \   / | \ /   |/
#    /|  \ /  | / \   |/
#    /|__/_\__|/___\__|/
#    /| \  /  | \ /   |/
#    /|_/__\__|_/_\___|/
#
function Malha5(nx::Int64=2,ny::Int64=2,etype=:truss2D;   
    Lx=4.0, Ly=4.0,  force=100E3, A=0.0 ,Ex=210E9,
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
           3 1; 
           3 2;
           4 1;
           4 2;
           6 1;
           6 2;
           7 1;
           7 2;
           9 1;
           9 2]

   # Generate the load information
   f = force
   nbc = [8 2 -f]

   # Todas as barras com seção circular e 2.5 cm de raio
   r = 2.5E-2
   A = pi*r^2

   # Material e geometria
   mat = [Material(Ex=Ex)]
   geom = [Geometry(A=A)]

   # Gera a malha e devolve um tipo Mesh2D
   Mesh2D(bmesh,mat,geom,hebc,nbc)

end


#
#  MALHA 6
#
#      ________________  
#    /| \   / | \ /   |
#    /|  \ /  | / \   |
#    /|__/_\__|/___\__| <---- P
#    /| \  /  | \ /   |
#    /|_/__\__|_/_\___|
#
function Malha6(nx::Int64=2,ny::Int64=2,etype=:truss2D;   
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
           4 1; 
           4 2;
           7 1;
           7 2]
             

   # Generate the load information
   f = force
   nbc = [6 1 -f]

   # Todas as barras com seção circular e 2.5 cm de raio
   r = 2.5E-2
   A = pi*r^2

   # Material e geometria
   mat = [Material(Ex=Ex)]
   geom = [Geometry(A=A)]

   # Gera a malha e devolve um tipo Mesh2D
   Mesh2D(bmesh,mat,geom,hebc,nbc)

end
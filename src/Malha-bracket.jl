
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

function Bracket2D(nx::Int64=3,ny::Int64=6,etype=:truss2D;   
    Lx=3.0, Ly=6.0,  force=10E3, A=0.0 ,Ex=210E9,
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
       25  1; 
       25  2]
              

# Generate the load information
f = force
nbc = [16 2 -f]

# Todas as barras com seção circular e 2.5 cm de raio
r = 2.5E-2
A = pi*r^2

# Material e geometria
mat = [Material(Ex=Ex)]
geom = [Geometry(A=A)]

# Gera a malha e devolve um tipo Mesh2D
Mesh2D(bmesh,mat,geom,ebc,nbc)

end

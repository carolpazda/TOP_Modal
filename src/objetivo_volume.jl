#
# Calcula o volume da estrutura
#
function Volume(m::Mesh,x::Vector)
   
    volume = 0.0
    
    # Loop for all elements
    for j=1:m.bmesh.ne
        
        # Area
        geo = m.geo_ele[j]
        A = m.geometries[geo].A
    
        # Length
        L = BMesh.Length(m.bmesh,j)  

        # Add to the total volume
        volume = volume + (L*A*x[j])
        
    end #j
    
    return volume
    
end

#
# Derivada parcial do volume em relação a x
#
function dVolume(m::Mesh,x::Vector)
    
    # Output 
    D = zeros(m.bmesh.ne)

    # Loop for all elements
    for j=1:m.bmesh.ne
    
        # Area
        geo = m.geo_ele[j]
        A = m.geometries[geo].A
        
        # Length
        L = BMesh.Length(m.bmesh,j)  

        # Local derivative
        D[j] = L*A
        
    end #j
    
    return D
    
end

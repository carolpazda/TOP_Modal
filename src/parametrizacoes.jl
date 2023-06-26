# Recebe a parametrização da rigidez e sua derivada
kparam(xe::Float64,p=penal)=xe^p
dkparam(xe::Float64,p=penal)=p*xe^(p-1)

# Recebe a parametrização da massa e sua derivada
mparam(xe::Float64,p=1.0)=xe^p
dmparam(xe::Float64,p=1.0)=p*xe^(p-1)

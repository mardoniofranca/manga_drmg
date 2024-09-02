ITensors.state(::StateName"Up", ::SiteType"S=3/2") = [1.0, 0.0, 0.0, 0.0]
ITensors.state(::StateName"Dn", ::SiteType"S=3/2") = [0.0, 0.0, 0.0, 1.0]
ITensors.state(::StateName"Up1", ::SiteType"S=3/2") = [0.0, 1.0, 0.0, 0.0]
ITensors.state(::StateName"Dn1", ::SiteType"S=3/2") = [0.0, 0.0, 1.0, 0.0]


function ITensors.op!(Op::ITensor,
                      ::OpName"Sz",
                      ::SiteType"S=3/2",
                      s::Index)
  Op[s'=>1,s=>1] = +3/2
  Op[s'=>2,s=>2] = +1/2
  Op[s'=>3,s=>3] = -1/2
  Op[s'=>4,s=>4] = -3/2
end

function ITensors.op!(Op::ITensor,
                      ::OpName"S+",
                      ::SiteType"S=3/2",
                      s::Index)
  Op[s'=>1,s=>2] = sqrt(3)
  Op[s'=>2,s=>3] = 2
  Op[s'=>3,s=>4] = sqrt(3)
end

function ITensors.op!(Op::ITensor,
                      ::OpName"S-",
                      ::SiteType"S=3/2",
                      s::Index)
  Op[s'=>2,s=>1] = sqrt(3)
  Op[s'=>3,s=>2] = 2
  Op[s'=>4,s=>3] = sqrt(3)
end

function ITensors.op!(Op::ITensor,
  ::OpName"Id",
  ::SiteType"S=3/2",
  s::Index)
Op[s'=>1,s=>1] = 1
Op[s'=>2,s=>2] = 1
Op[s'=>3,s=>3] = 1
Op[s'=>4,s=>4] = 1
end


#################### quantum numbers ##########

function ITensors.space(::SiteType"S=3/2";
  conserve_qns=false)
  if conserve_qns
    return [QN("Sz",3)=>1,QN("Sz",1)=>1, QN("Sz",-1)=>1,QN("Sz",-3)=>1]
  end

return 4
end



################## Setting Hamiltonian ###############

function Set_Hamiltonian(N,J_1,J_2,D,E,h0,sites)::MPO
  ampo = OpSum()
  for j=1:N-1
    ampo += 0.5*J_1,"S+",j,"S-",j+1
    ampo += 0.5*J_1,"S-",j,"S+",j+1
    ampo += J_1,"Sz",j,"Sz",j+1
  end

  for j=1:N-2
    ampo += 0.5*J_2,"S+",j,"S-",j+2
    ampo += 0.5*J_2,"S-",j,"S+",j+2
    ampo += J_2,"Sz",j,"Sz",j+2  
  end

  for j=1:N
    ampo += D,"Sz",j,"Sz",j

    # Sx
    ampo += 0.5*E,"S+",j
    ampo += 0.5*E,"S-",j
    #Sy
    ampo += -0.5*1im*E,"S+",j
    ampo += +0.5*1im*E,"S+",j

    ######### field

    ampo += h0*j,"Sz",j

  end
  
  HH=MPO(ampo,sites)
  return HH
end


######################## Imbalance #####################

function Imbalance(M1,M2,N)
  S=0
  for i=1:N 
    S+=M1[i,i]*M2[i,i]
  end
  return real(4*S/(9*N))
end


##################### ENTROPY #########################

function VN_entropy(psi::MPS, b::Int)
  s = siteinds(psi)  
  orthogonalize!(psi, b)
  _,S = svd(psi[b], (linkind(psi, b-1), s[b]))
  SvN = 0.0
  for n in 1:dim(S, 1)
    p = S[n,n]^2
    SvN -= p * log(p)
  end
  return SvN
end


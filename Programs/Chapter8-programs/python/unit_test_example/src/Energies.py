import PairingHamiltonian as ham
import PairingBasisHardCoded as basis
# import PairingBasisGenerated as pbg
# basis = pbg.PairingBasisGen(8,4)

def eps(holes, particles):
    E = 0
    for h in holes:
        p, s = basis.states[h]
        E += (p-1)
    for p in particles:
        p, s = basis.states[p]
        E -= (p-1)
    return E


# Diagram 1
def d1():
    sum = 0
    for a in basis.above_fermi:
        for b in basis.above_fermi:
            for i in basis.below_fermi:
                for j in basis.below_fermi:
                    sum += 0.25*ham.assym(a,b,i,j)*ham.assym(i,j,a,b)/eps((i,j),(a,b))
    return sum



# Diagram 3
def d3():
    sum = 0
    for a in basis.above_fermi:
       for b in basis.above_fermi:
           for c in basis.above_fermi:
               for i in basis.below_fermi:
                   for j in basis.below_fermi:
                       for k in basis.below_fermi:
                           sum += ham.assym(i,j,a,b)*ham.assym(a,c,j,k)*ham.assym(b,k,c,i)/eps((i,j),(a,b))/eps((k,j),(a,c))
    return sum

# Diagram 4
def d4():
    sum = 0
    for a in basis.above_fermi:
        for b in basis.above_fermi:
            for c in basis.above_fermi:
                for d in basis.above_fermi:
                    for i in basis.below_fermi:
                        for j in basis.below_fermi:
                            sum += 0.125*ham.assym(i,j,a,b)*ham.assym(a,b,c,d)*ham.assym(c,d,i,j)/eps((i,j),(a,b))/eps((i,j),(c,d))
    return sum

# Diagram 5
def d5():
    sum = 0
    for a in basis.above_fermi:
        for b in basis.above_fermi:
            for i in basis.below_fermi:
                for j in basis.below_fermi:
                    for k in basis.below_fermi:
                        for l in basis.below_fermi:
                            sum += 0.125*ham.assym(i,j,a,b)*ham.assym(k,l,i,j)*ham.assym(a,b,k,l)/eps((i,j),(a,b))/eps((k,l),(a,b))
    return sum

# Diagram 8 
def d8():
    sum = 0
    for a in basis.above_fermi:
        for b in basis.above_fermi:
            for i in basis.below_fermi:
                for j in basis.below_fermi:
                    for k in basis.below_fermi:
                        sum -= 0.5*ham.assym(i,j,a,b)*ham.assym(a,b,i,k)*ham.f(k,j)/eps((i,j),(a,b))/eps((i,k),(a,b))
    return sum

# Diagram 9 
def d9():
    sum = 0
    for a in basis.above_fermi:
        for b in basis.above_fermi:
            for c in basis.above_fermi:
                for i in basis.below_fermi:
                    for j in basis.below_fermi:
                        sum += 0.5*ham.assym(i,j,a,b)*ham.assym(a,c,i,j)*ham.f(b,c)/eps((i,j),(a,b))/eps((i,j),(a,c))
    return sum

def corr_mbpt2():
    sum = d1()
    return sum

def corr_mbpt3():
    sum = d1() + d3() + d4() + d5()
    return sum

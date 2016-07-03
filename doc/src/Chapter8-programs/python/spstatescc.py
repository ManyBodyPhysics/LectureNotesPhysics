from numpy import *

nmax =2
nshell = 3*nmax*nmax
count = 1
tzmin = 1

print "Symmetric nuclear matter:"  
print "a, nx,   ny,   nz,   sz,   tz,   nx^2 + ny^2 + nz^2"
for n in range(nshell): 
    for nx in range(-nmax,nmax+1):
         for ny in range(-nmax,nmax+1):
            for nz in range(-nmax, nmax+1):  
                for sz in range(-1,1+1):
                    tz = 1
                    for tz in range(-tzmin,tzmin+1):
                        e = nx*nx + ny*ny + nz*nz
                        if e == n:
                            if sz != 0: 
                                if tz != 0: 
                                    print count, "  ",nx,"  ",ny, "  ",nz,"  ",sz,"  ",tz,"         ",e
                                    count += 1
                                    
                                    
nmax =1
nshell = 3*nmax*nmax
count = 1
tzmin = 1
print "------------------------------------"
print "Neutron matter:"                                    
print "a, nx,   ny,   nz,   sz,    nx^2 + ny^2 + nz^2"
for n in range(nshell): 
    for nx in range(-nmax,nmax+1):
         for ny in range(-nmax,nmax+1):
            for nz in range(-nmax, nmax+1):  
                for sz in range(-1,1+1):
                    e = nx*nx + ny*ny + nz*nz
                    if e == n:
                        if sz != 0: 
                            print count, "  ",nx,"  ",ny, "  ",sz,"  ",tz,"         ",e
                            count += 1     

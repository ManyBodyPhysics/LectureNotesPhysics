import sys
sys.path.append('../src/')

import PairingBasisGenerated as pbg
import PairingBasisHardCoded as pbhc

basis = pbg.PairingBasisGen(8,4)

flag = 0
for elem in range(len(pbhc.states)):
    if( pbhc.states[elem] != basis.states[elem] ):
        print pbhc.states[elem], basis.states[elem], "aren't equal!"
        flag = 1

if( flag == 0 ):
    print "test passed!"
else:
    print "test failed."

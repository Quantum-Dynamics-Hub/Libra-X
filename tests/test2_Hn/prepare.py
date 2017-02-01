
a = 0.708  # equilibrium H-H distance in H2 molecule (just bonded terms of UFF)
N = 10

for n in xrange(N):
    print "H 1.0     0    0    %8.5f " % (n*a)

for n in xrange(N):

    if n==0:
        print "CONECT %5i %5i" % (n+1, n+1+1)
    elif n==N-1:
        print "CONECT %5i %5i" % (n+1, n-1+1)
    else:
        print "CONECT %5i %5i %5i" % (n+1, n-1+1, n+1+1)



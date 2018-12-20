import matplotlib.pyplot as plt
import subprocess, os
import numpy as np

def bench( distribution, MDPC ):
    ns = [ 1e5, 3e5, 1e6, 3e6, 1e7, 3e7 ] # 

    # internal
    ts = []
    tm = []
    for n in ns:
        result = subprocess.check_output( "nsmake run -DMDPC=pouet benchmarks/bench_traversal_2D.cpp {} {} {}".format( distribution, MDPC, int( n ) ).split() )
        spl = str( result ).replace( "\\n", " " ).split()
        add = spl[ spl.index( "[heap]" ) + 1 ].split( "-" )
        mem = int( add[ 1 ], 16 ) - int( add[ 0 ], 16 )
        ts.append( float( spl[ -2  ] ) )
        tm.append( mem )
        print( ts, tm )

    # CGAL
    cs = []
    cm = []
    for n in ns:
        result = subprocess.check_output( "./scalability_2D_cgal.exe {} {}".format( distribution, int( n ) ).split() )
        spl = str( result ).replace( "\\n", " " ).split()
        add = spl[ spl.index( "[heap]" ) + 1 ].split( "-" )
        mem = int( add[ 1 ], 16 ) - int( add[ 0 ], 16 )
        cs.append( float( spl[ -2  ] ) )
        cm.append( mem )
        print( cs, cm )

    # print( ts, cs )
    print( "time:", cs[ -1 ] / ts[ -1 ] )
    print( "memo:", cm[ -1 ] / tm[ -1 ] )
    print( "times:", np.array( cs ) / np.array( ts ) )

    plt.loglog( ns, ts )
    plt.loglog( ns, cs )
    plt.legend( [ "internal", "cgal" ] )
    plt.savefig( "bench_traversal_speed_{}.png".format( distribution ) )
    # plt.show()
    plt.clf()

    plt.loglog( ns, tm )
    plt.loglog( ns, cm )
    plt.legend( [ "internal (mem)", "cgal (mem)" ] )
    plt.savefig( "bench_traversal_mem_{}.png".format( distribution ) )
    # plt.show()
    plt.clf()

os.system( "g++ -c -o scalability_2D_cgal.o -march=native -O3 -lCGAL benchmarks/bench_traversal_2D_CGAL.cpp" )
os.system( "g++ -o scalability_2D_cgal.exe scalability_2D_cgal.o -lCGAL -lgmp" )

# bench( "lines"  , 38 )
bench( "regular", 11 )
# bench( "random" , 18 )


import matplotlib.pyplot as plt
import numpy as np

txt = """
random 10000 True (18, 1) 0.037519  0.004642
random 100000 True (14, 1) 0.384125 0.047407
random 1000000 True (10, 1) 3.87565 0.6074
random 10000000 True (20, 1) 39.4562 7.24189

regular 10000 True (18, 1) 0.0376 0.004007
regular 100000 True (10, 1) 0.378673 0.038261
regular 1000000 True (10, 1) 3.86604 0.36494
regular 10000000 True (20, 1) 39.2951 4.20647

random 10000 False (30, 150) 0.037954 0.012483
random 100000 False (50, 96) 0.391179 0.135544
random 1000000 False (30, 150) 3.9446 2.01066

regular 10000 False (12, 1) 0.037168 0.004141
regular 100000 False (12, 0.5) 0.379598 0.038869
regular 1000000 False (12, 2) 3.87423 0.385696

split 10000 True (10, 20) 0.037883 0.004652
split 100000 True (10, 20) 0.384111 0.052311
split 1000000 True (20, 20) 3.87505 0.634518

split 10000 False (35, 460) 0.037543 0.029782

"""

res = {}
for l in txt.split("\n"):
    l = l.replace( "(", "" ).replace( ")", "" ).replace( "  ", " " )
    s = l.split( " " )
    if len( s ) < 5:
        continue
    if s[ 2 ] == "True":
        t = s[ 0 ] + '_vor'
    else:
        t = s[ 0 ] + '_pow'
    if not ( t in res ):
        res[ t ] = {
            "size": [],
            "cgal": [],
            "pdgm": []
        }
    print( s )
    res[ t ][ "size" ].append( int  ( s[ 1 ] ) )
    res[ t ][ "cgal" ].append( float( s[ 5 ] ) )
    res[ t ][ "pdgm" ].append( float( s[ 6 ] ) )

for t, v in res.items():
    print( t, v )
    plt.clf()
    plt.loglog( v[ "size" ], v[ "cgal" ] )
    plt.loglog( v[ "size" ], v[ "pdgm" ] )
    plt.legend( [ "cgal", "pdgm" ] )
    plt.xlabel( "nb diracs" )
    plt.ylabel( "time (s)" )
    plt.grid((True,True))
    plt.title( t )
    # plt.show()
    plt.savefig( "vtk/" + t + ".png" )

    print( t, np.array( v[ "cgal" ] ) / np.array( v[ "pdgm" ] ) )

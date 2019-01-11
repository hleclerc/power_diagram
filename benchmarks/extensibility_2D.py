import os, sys, subprocess

# nws = [ int( 1e4 ), int( 1e5 ), int( 1e6 ) ]
# distributions = [ "random", "regular", "split" ]
voronois = [ False, True ]

distributions = [ "split" ]
voronois = [ False ]
nws = [ int( 1e4 ), int( 1e5 ), int( 1e6 ) ]

# get optimal weights for each distribution
for distribution in distributions:
    for nw in nws:
        for voronoi in voronois:
            if voronoi:
                inp = distribution + " --nb-diracs=" + str( nw )
            else:
                inp = "vtk/" + distribution + "_" + str( nw ) + ".xyw"
                if not os.path.exists( inp ):
                    os.system( "nsmake run -g3 benchmarks/bench_solve_2D.cpp --distribution=" + distribution + " --nb-diracs=" + str( nw ) + " --max-dirac-per-cell=35 --max-delta-weight=1e40 --output=" + inp )

            timing_pdgm = 1e6
            timing_cgal = 1e6
            best_conf = ()
            for max_dirac_per_cell in [ 35 ]: # 5, 10, 30, 50
                lrmdw = [ 460, 1000 ] # 0.5, 1, 2, 4, 8, 16
                for rmdw in lrmdw:
                    if voronoi and rmdw != lrmdw[ 0 ]:
                        continue
                    mdw = rmdw / nw
                    cmd = "nsmake run -g3 benchmarks/bench_solve_2D.cpp --max-iter=1 --distribution=" + inp + " --max-dirac-per-cell=" + str( max_dirac_per_cell ) + " --max-delta-weight=" + str( mdw )
                    res = subprocess.run( cmd.split( " " ), stdout = subprocess.PIPE )
                    for l in res.stdout.decode( "utf8" ).split( "\n" ):
                        if l.startswith( "timings_der" ):
                            p = float( l.split( " " )[ -1 ] )
                            if timing_pdgm > p:
                                best_conf = ( max_dirac_per_cell, rmdw )
                                timing_pdgm = p
                        if l.startswith( "timings_cgal" ):
                            timing_cgal = min( timing_cgal, float( l.split( " " )[ -1 ] ) )
                                

            print( distribution, nw, voronoi, best_conf, timing_cgal, timing_pdgm )


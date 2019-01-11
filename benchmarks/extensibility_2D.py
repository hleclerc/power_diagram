import os, sys, subprocess

distributions = [ "random" ]
nws = [ int( 1e4 ), int( 1e5 ) ]
# , int( 1e6 )

# get optimal weights for each distribution
# for distribution in distributions:
#     for nw in nws:
#         mdw = 60 / nw
#         os.system( "nsmake run -g3 benchmarks/bench_solve_2D.cpp --distribution " + distribution + " --nb-diracs " + str( nw ) + " --max-dirac-per-cell 35 --max-delta-weight " + str( mdw ) + " --output vtk/" + distribution + "_" + str( nw ) + ".xyw" )

# get timings
for distribution in distributions:
    for nw in nws:
        for voronoi in [ True, False ]:
            if voronoi:
                input = "--distribution " + distribution + " --nb-diracs " + str( nw )
            else:
                input = "--distribution vtk/" + distribution + "_" + str( nw ) + ".xyw"

            timing_pdgm = 1e6
            timing_cgal = 1e6
            for max_dirac_per_cell in [ 10, 20, 30, 40, 50, 60 ]:
                for mdw in [ 6 / nw, 12 / nw, 24 / nw, 60 / nw, 100 / nw, 200 / nw ]:
                    cmd = "nsmake run -g3 benchmarks/bench_solve_2D.cpp --max-iter 1 " + input + " --max-dirac-per-cell " + str( max_dirac_per_cell ) + " --max-delta-weight " + str( mdw )
                    res = subprocess.run( cmd.split( " " ), stdout = subprocess.PIPE )
                    for l in res.stdout.decode( "utf8" ).split( "\n" ):
                        if l.startswith( "timings_der" ):
                            timing_pdgm = min( timing_pdgm, float( l.split( " " )[ -1 ] ) )
                        if l.startswith( "timings_cgal" ):
                            timing_cgal = min( timing_cgal, float( l.split( " " )[ -1 ] ) )
                    print( timing_pdgm )

            print( distribution, nw, voronoi, timing_cgal, timing_pdgm )


# Real Data
Typical config files for use with real data for unversions.


run-escript ~/MINTpy/bin/getRealDataEDI.py

run-escript ~/MINTpy/bin/mesh/mkMeshEDI.py configFirstOrder

run-escript ~/MINTpy/bin/mesh/mkFly mesh/shallow

run-escript -t10 ~/MINTpy/bin/MT_1stOrderScaled.py configFirstOrder

run-escript -t10 ~/MINTpy/bin/MT_1stOrderScaled.py configSecondOrder

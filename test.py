Helium = Abel_Noble_Gas(5190.0,2.83e-3,4e-3,5./3,diameter=0.5e-3,P0=1.0e7)
Helium.SolveNozzleProperty()
Helium.Birch87Solver()
Helium.MachDiskSolver()
Helium.SlippingRegionSolverT()
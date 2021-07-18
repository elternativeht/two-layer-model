from Twolayer_model import Abel_Noble_Gas

Helium = Abel_Noble_Gas(cp=5190.0, b = 2.83e-3, M=4e-3, gamma=5./3, nozzle_diameter=0.5e-3, P0=1.0e7)
Helium.SolveNozzleProperty()
Helium.Birch87Solver()
Helium.MachDiskSolver()
Helium.SlippingRegionSolverT()
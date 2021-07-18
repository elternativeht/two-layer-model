from math import pi
from scipy.optimize import fsolve
import numpy as np
class Abel_Noble_Gas(object):
    def __init__(self,cp,b,M,gamma,nozzle_diameter=0.5e-3,P0=1.5e7,T0=300.0):
        R = 8.314
        #M in kg/mole
        assert(b>0 and M>0 and P0>0 and T0>0)
        self.cp = cp
        self.gamma = gamma
        self.Rg = R/M
        self.b = b
        self.p0 = P0
        self.T0 = T0
        self.rho_0 = self.Rg * self.T0 / self.p0 + self.b
        self.rho_0 = 1.0/self.rho_0
        self.diameter = nozzle_diameter
        self.NozzleArea = 0.25*pi*(self.diameter**2)
        self.MassFlow = 0.0
    def PropertyOutput(self):
        print("Current pressure: "+str(self.p))
        print("Current temperature: "+str(self.T))
        print("Current density: "+str(self.rho))
    def SolveNozzleProperty(self,isentropic=True):
        p0 = self.p0
        Rg = self.Rg
        T0 = self.T0
        rho_0 = self.rho_0
        cp = self.cp
        b = self.b
        gamma = self.gamma
        def Eqs(proplist):
            nonlocal T0, rho_0,b,cp,gamma,Rg,p0
            p1,T1,rho_1,U1 = proplist
            return [
                cp*T0+b*p0-cp*T1-b*p1-(U1**2)/2,
                U1-(1/(1-b*rho_1))*((gamma*Rg*T1)**0.5),
                p1*(1-b*rho_1)/rho_1-Rg*T1,
                p1*((1-b*rho_1)/rho_1)**(gamma) - p0*((1.0-b*rho_0)/rho_0)**gamma
            ]
        ResultList = fsolve(Eqs,[0.5*self.p0,0.8*self.T0,0.6*self.rho_0,250.0])
        self.p1,self.T1,self.rho_1,self.U1 = ResultList
        self.MassFlow = self.U1*self.rho_1*self.NozzleArea
        print('''
        Choked nozzle Properties:
        -----------------------------------
        Stagnation Pressure:   %f MPa
        Nozzle Diameter:       %f mm
        Nozzle Temperature:    %f K
        Nozzle Velocity:       %f m/s
        Nozzle Density:        %f kg/m^3
        Nozzle mass flow rate: %f kg/s
        -----------------------------------
        ''' %(self.p0/1.0e6,self.diameter*1.0e3,self.T1,self.U1,self.rho_1,self.MassFlow))
    def Birch87Solver(self,p_intfy=101325.,T_B2=300.0):
        def Eqs_B87(param_list):
            nonlocal p_intfy
            rho_B2,U_B2,A_B2 = param_list
            v_B2 = 1./rho_B2
            return [
                p_intfy*(v_B2-self.b)-self.Rg*T_B2,
                (self.p1-p_intfy)*self.NozzleArea+self.MassFlow*self.U1-A_B2*(U_B2**2)/v_B2,
                U_B2*A_B2/v_B2-self.MassFlow
            ]
        he_init_value = [0.50,500.0,1e-4]
        h2_init_value = [8.11977575e-02,2.11509767e+03,1.02122048e-05]
        h2_init_value_new = [8.30938830e-02,2.06612138e+03,1.11813534e-03]
        resultlist = fsolve(Eqs_B87,h2_init_value_new)
        print(resultlist)
        self.rho_B2,self.U_B2,self.A_B2 = resultlist
        self.Birch87_Diameter = (self.A_B2*4/pi)**0.5
        print('''
        Birch87 Boundary Conditions:
        -----------------------------
        Notional diameter: %f mm
        Notional Velocity: %f m/s
        ''' % (self.Birch87_Diameter*1000,self.U_B2))
    def MachDiskSolver(self,p_infty=101325.0,gas_property='he'):
        self.p2a = 0.0
        self.T2a = 0.0
        self.rho_2a = 0.0
        self.U2a = 0.0
        self.p2b = p_infty
        self.T2b = 0.0
        self.rho_2b = 0.0
        self.U2b = 0.0
        ratio = 0.18
        self.MachDiskDiameter = 0.35*((self.p0/p_infty)**0.5)*self.diameter
        self.MachDiskCutArea = 0.25*pi*(self.MachDiskDiameter**2)
        self.MachDiskLocation = 0.67*((self.p0/p_infty)**0.5)*self.diameter
        tempTC = 1.0-(self.gamma+1)*(((self.gamma+1)/(self.gamma-1))**(-0.5))/self.gamma
        self.MachDiskDiameter_T = 0.954*self.MachDiskLocation*tempTC**0.5
        #print(self.MachDiskLocation)
        #print(tempTC)
        #print(self.MachDiskDiameter_T)
        self.MachDiskCutArea_T = 0.25*pi*(self.MachDiskDiameter_T**2)
        if self.p0>5.0e6:
            pressurelist = [1.0e6,3.0e6,5.0e6,7.0e6,self.p0]
        else:
            pressurelist = [5e5,7e5,9e5,1e6,2.0e6,self.p0]
        he_init_value_ = [1.70324882e+03,1.77509048e+01,4.61584098e-02,1.71495072e+03,2.80962928e+02,1.73422070e-01,4.56455155e+02]
        h2_init_value_ = [5.07478572e+03,6.63019520e+01,1.84098341e-02,2.58672039e+03,2.89287039e+02,8.42027490e-02,5.65552713e+02]
        _init_value_ = he_init_value_[:] if gas_property == 'he' else h2_init_value_
        #_init_value_ = [p_infty*0.5*ratio,self.T0*ratio,0.1*self.rho_0*ratio,130.0*ratio,self.T0*ratio,0.1*self.rho_0*ratio,130.0*ratio]
        for cur_p0 in pressurelist:
            cur_T0 = self.T0
            cur_v0 = self.Rg*cur_T0/cur_p0+self.b
            cur_rho0 = 1./cur_v0
            def Eq_mach(proplist_):
                nonlocal cur_p0,cur_T0,cur_rho0
                p2a,T2a,rho_2a,U2a,T2b,rho_2b,U2b = proplist_
                v2a =1./rho_2a
                v2b = 1./rho_2b
                v0 = 1./cur_rho0
                return [
                    p2a*(v2a-self.b)-self.Rg*T2a,
                    self.cp*cur_T0+self.b*cur_p0-self.cp*T2a-self.b*p2a-(U2a**2)/2,
                    p2a*(v2a-self.b)**(self.gamma)-cur_p0*(v0-self.b)**(self.gamma),
                    U2a/v2a-U2b/v2b,
                    self.p2b*(v2b-self.b)-self.Rg*T2b,
                    p2a+(U2a**2)/v2a-self.p2b-(U2b**2)/v2b,
                    self.cp*T2b+self.b*self.p2b+(U2b**2)/2.-self.cp*cur_T0-self.b*cur_p0   
                ]
            #init_value_list = [1.70324882e+03,1.77509048e+01,4.61584098e-02,1.71495072e+03,2.80962928e+02,1.73422070e-01,4.56455155e+02]
            cur_ans = fsolve(Eq_mach,_init_value_)
            if cur_p0 != self.p0:
                _init_value_ = cur_ans
            else:
                print(cur_ans)
                self.p2a,self.T2a,self.rho_2a,self.U2a,self.T2b,self.rho_2b,self.U2b = cur_ans
                self.MachDiskMassFlow = self.rho_2b*self.U2b*self.MachDiskCutArea
                self.MachDiskMassFlow_T = self.rho_2b*self.U2b*self.MachDiskCutArea_T
                v2b = 1./self.rho_2b
                Vsound = v2b/(v2b-self.b)*(self.gamma*self.Rg*self.T2b)**0.5
                self.MachNumber = self.U2b/Vsound
                #print(cur_ans)
        print('''
        Post Mach Disk Properties:
        ------------------------------
        Post Mach Disk Temperature: %f K.
        Post Mach Disk Density:     %f kg/m^3
        Post Mach Disk Velocity:    %f m/s
        Post Mach Disk Diameter:    %f mm
        Post Mach Disk mass flowrate: %f kg/s
        Post Mach Disk Mach number: %f 
        Mach Disk Location:         %f mm
        ------------------------------
        Post Mach Disk Diameter(T): %f mm
        Post Mach Disk massflow(T): %f   kg/s
        ''' % (self.T2b,self.rho_2b,self.U2b,self.MachDiskDiameter*1000,self.MachDiskMassFlow,self.MachNumber,self.MachDiskLocation*1000,self.MachDiskDiameter_T*1000,self.MachDiskMassFlow_T))
    def SlippingRegionSolver(self,p_infty=101325.):
        self.SlippingRegionThickness = 0.30*self.diameter*(self.p0/p_infty)**0.5
        self.SlippingRegionArea = pi*(self.MachDiskDiameter/2.+self.SlippingRegionThickness)**2-self.MachDiskCutArea
        self.p3 = 101325.
        cp_air = 1005.
        T_infty_air = 300.
        he_init_guess =[240.,0.100,1.2,900.,40000.,101325-40000.]
        init_guess = he_init_guess[:]
        T_Energy = (self.MassFlow-self.MachDiskMassFlow)*(self.cp*self.T0+self.b*self.p0)
        def Eqs_SR(param_list):
            nonlocal p_infty,cp_air,T_infty_air,T_Energy
            Rg_air = 8314/29.0
            T3,rho_he3,rho_air3,U3,p_he3,p_air3 = param_list
            v_he3 = 1./rho_he3
            v_air3 = 1./rho_air3
            return [
                p_he3 + p_air3 - self.p3,
                p_he3*(v_he3-self.b)-self.Rg*T3,
                p_air3*v_air3-Rg_air*T3,
                U3*self.SlippingRegionArea/v_he3+self.MachDiskMassFlow-self.MassFlow,
                (self.p0-p_infty)*self.NozzleArea+self.MassFlow*self.U1-self.MachDiskMassFlow*self.U2b-(U3**2)*self.SlippingRegionArea*(1/v_he3+1/v_air3),
                U3*self.SlippingRegionArea*(self.cp*T3+self.b*p_infty+U3**2/2.)/v_he3+U3*self.SlippingRegionArea*(cp_air*(T3-T_infty_air)+U3**2/2.)/v_air3-T_Energy
            ]
        solve_result = fsolve(Eqs_SR,init_guess)
        self.T3,self.rho_he3,self.rho_air3,self.U3,self.p_he3,self.p_air3 = fsolve(Eqs_SR,init_guess)
        self.SlippingRegionHeMassFlow = self.rho_he3*self.U3*self.SlippingRegionArea
        self.SlippingRegionAirMassFlow = self.rho_air3*self.U3*self.SlippingRegionArea
        self.SlippingRegionTotalMassFlow = self.SlippingRegionHeMassFlow + self.SlippingRegionAirMassFlow
        self.SlippingRegionHeMassFraction = self.SlippingRegionHeMassFlow/self.SlippingRegionTotalMassFlow
        print('''
        Mixing Layer Properties:
        --------------------------------
        Thickness:      %f mm
        Pressure:       %f Pa
        Temperature:    %f K
        Velocity:       %f m/s
        Mass flow rate: %f kg/s
        He Mass Flow:   %f kg/s
        He Fraction:    %f
        ''' %(self.SlippingRegionThickness*1000,p_infty,self.T3,self.U3,self.SlippingRegionTotalMassFlow,self.SlippingRegionHeMassFlow,self.SlippingRegionHeMassFraction))
    def SlippingRegionSolverT(self,p_infty=101325.):
        temp_cc = (self.p0/p_infty)**((self.gamma-1)/self.gamma)*(1+(self.gamma-1)/2)
        self.SlippingRegionThicknessT = 0.135*self.MachDiskLocation*(1+1./temp_cc)
        self.SlippingRegionAreaT = pi*(self.MachDiskDiameter_T/2.+self.SlippingRegionThicknessT)**2-self.MachDiskCutArea_T
        self.p3 = 101325.
        cp_air = 1005.
        T_infty_air = 300.
        he_init_guess =[240.,0.100,1.2,900.,40000.,101325-40000.]
        new_he_init_guess = [9.43562665e+01,4.65699850e-01,3.64927782e-01,1.19400673e+03,9.14533510e+04,9.87164895e+03]
        new_h2_init_guess =[1.35619752e+02,1.68657692e-01,1.57330526e-01,1.65748130e+03,9.52078664e+04,6.11713356e+03]
        init_guess = new_h2_init_guess[:]
        T_Energy = (self.MassFlow-self.MachDiskMassFlow_T)*(self.cp*self.T0+self.b*self.p0)
        def Eqs_SR(param_list):
            nonlocal p_infty,cp_air,T_infty_air,T_Energy
            Rg_air = 8314/29.0
            density_air_out = p_infty/(Rg_air*T_infty_air)
            T3,rho_he3,rho_air3,U3,p_he3,p_air3 = param_list
            #mlayerarea = pi*(self.MachDiskDiameter_T/2.+mlayer_thickness)**2-self.MachDiskCutArea_T
            v_he3 = 1./rho_he3
            v_air3 = 1./rho_air3
            return [
                p_he3 + p_air3 - self.p3,
                p_he3*(v_he3-self.b)-self.Rg*T3,
                p_air3*v_air3-Rg_air*T3,
                U3*self.SlippingRegionAreaT/v_he3+self.MachDiskMassFlow_T-self.MassFlow,
                (self.p0-p_infty)*self.NozzleArea+self.MassFlow*self.U1-self.MachDiskMassFlow_T*self.U2b-(U3**2)*self.SlippingRegionAreaT*(1/v_he3+1/v_air3),
                U3*self.SlippingRegionAreaT*(self.cp*T3+self.b*p_infty+U3**2/2.)/v_he3+U3*self.SlippingRegionAreaT*(cp_air*(T3-T_infty_air)+U3**2/2.)/v_air3-T_Energy
                #mlayer_thickness/self.MachDiskLocation-0.135*(1+density_air_out/())
            ]
        q = fsolve(Eqs_SR,init_guess)
        #print(q)
        self.T3,self.rho_he3,self.rho_air3,self.U3,self.p_he3,self.p_air3 = q
        self.SlippingRegionHeMassFlowT = self.rho_he3*self.U3*self.SlippingRegionAreaT
        #print(self.MassFlow)
        #print(self.MachDiskMassFlow_T)
        #print(self.SlippingRegionHeMassFlowT)
        self.SlippingRegionAirMassFlowT = self.rho_air3*self.U3*self.SlippingRegionAreaT
        self.SlippingRegionTotalMassFlowT = self.SlippingRegionHeMassFlowT + self.SlippingRegionAirMassFlowT
        self.SlippingRegionHeMassFractionT = self.SlippingRegionHeMassFlowT/self.SlippingRegionTotalMassFlowT
        print('''
        Mixing Layer Properties:
        --------------------------------
        Thickness:      %f mm
        Pressure:       %f Pa
        Temperature:    %f K
        Velocity:       %f m/s
        Mass flow rate: %f kg/s
        He Mass Flow:   %f kg/s
        He Fraction:    %f
        ''' %(self.SlippingRegionThicknessT*1000,p_infty,self.T3,self.U3,self.SlippingRegionTotalMassFlowT,self.SlippingRegionHeMassFlowT,self.SlippingRegionHeMassFractionT))
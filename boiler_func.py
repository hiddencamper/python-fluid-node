import iapws
import iapws._iapws97Constants
import numpy as np
import matplotlib.pyplot as plt
import math

class Cylindrical_Boiler:

    gasConstant = 8.3145
    molarMassOfSteam = 18.01528 / 1000 # kg/ mol
    PSI_TO_PASCAL = 6895
    gravity = 9.81
    carbonSteelCp = .49 # jk / (kg * K)

    def idealGas_Mass(self):
        returnMass = (self.PressurePa * self.SteamVol * self.molarMassOfSteam)/(self.TempK * self.gasConstant)
        return returnMass

    def idealGas_PressurePa(self):
        self.PressurePa = (self.Steam_Mass/self.molarMassOfSteam) * self.gasConstant * self.TempK / self.SteamVol #Pa
        self.PressureMPa = self.PressurePa / 1e6 #MPa
        self.PressPSIG = self.PressurePa / self.PSI_TO_PASCAL #PSIG
        return self.PressurePa

    def Temp_Celsius(self):
        return self.TempK - 273.15

    def __init__(self, initradius, initheight, initH20Vol, initTempK, initPress, initMass, initdeltaT):
        self.radius = initradius
        self.height = initheight
        self.Area = math.pi*self.radius*self.radius
        self.Volume = self.Area*self.height
        self.WaterVol = initH20Vol
        self.TempK = initTempK
        self.PressureMPa = initPress
        self.InternalMass = initMass
        state = iapws.IAPWS97(T=self.TempK, P=self.PressureMPa)
        self.Enthalpy = state.h # initial enthalpy of system is in equilibrium
        self.Liquid_Mass = self.WaterVol * state.Liquid.rho
        if state.region == 4:
            self.DropletMass = self.WaterVol * state.Vapor.rho #droplets if we are boiling/2 phase
            self.DropletVol = self.DropletMass * state.Vapor.v
        else:
            self.DropletMass = 0 
            self.DropletVol = 0
        self.LiquidSpaceVol = self.WaterVol + self.DropletVol # Volume of the liquid space, water + droplet volume
        self.SteamVol = self.Volume - self.LiquidSpaceVol #steam volume is the remainder of volume
        self.PressurePa = self.PressureMPa / 1e6
        self.PressPSIG = self.PressurePa / 6895
        self.Steam_Mass = self.idealGas_Mass()
        self.Total_Mass_LiquidSpace = self.Liquid_Mass + self.DropletMass
        self.StepSize = 1
        self.heatingPower = 0
        self.FeedwaterFlow = 0
        self.SteamFlow = 0
        self.LiquidHeight = self.LiquidSpaceVol / self.Area
        self.InternalMassEnergy = self.InternalMass * self.carbonSteelCp * self.TempK
        self.runFW = True
        self.runSteamflow = True
        self.Quality = state.x
        self.LastDropletFlow = 0
        self.LastDropletVel = 0
        self.deltaT = initdeltaT
        self.LastSteamDemand = 0
        self.LastSwellHeight = 0
        self.LastLiquidVol = self.LiquidSpaceVol

    def Set_Heating_Power(self, newPower):
        self.heatingPower = newPower

    def Droplet_Travel(self,state):        
        if state.region != 4:
            return
#        Pbot = state.rho * self.gravity * self.LiquidHeight + self.PressurePa
#        PbotPSIG=Pbot/self.PSI_TO_PASCAL
#        top = abs(Pbot - self.PressurePa) * self.Area * self.Area
#        bot = 8 * math.pi * state.Liquid.mu * self.LiquidHeight
#        dropletFlow = top/bot

        Pbot = state.rho * self.gravity * self.LiquidHeight
        dropVel = (2 * Pbot / state.rho) #TODO: /10 to help slow down process
        #TODO: Need to come up with a new way to model droplet velocity
        #TODO: Droplet flow in a unit of time should be velocity * dropletmass. 
        #if we make x mass of droplets per second, and it takes y seconds for them all to 
        #rise up to the liquid level, then that should set the droplet residence time
        #for the swell effect to work, we need the droplets in the water to expand
        #due to lowering pressure.

        dropletFlow = dropVel * self.DropletMass / self.LiquidSpaceVol
        
        dropletFlow = self.DropletMass #TODO:Comment this out

        if dropletFlow > self.DropletMass:
            dropletFlow = self.DropletMass
            print("droplet_flow_max")
        self.DropletMass -= dropletFlow
        self.Total_Mass_LiquidSpace -= dropletFlow
        self.Enthalpy = (self.Liquid_Mass + self.DropletMass) * state.h/self.Total_Mass_LiquidSpace 
        self.Steam_Mass += dropletFlow
        self.LastDropletFlow = dropletFlow
        self.LastDropletVel = dropVel

    def injectFW(self,FW_mass_flow_ks_sec, FW_enthalpy):
        self.RunFW = True
        self.Liquid_Mass += FW_mass_flow_ks_sec
        self.Enthalpy = (self.Total_Mass_LiquidSpace * self.Enthalpy + FW_mass_flow_ks_sec*FW_enthalpy)/(self.Total_Mass_LiquidSpace+FW_mass_flow_ks_sec) 
        self.Total_Mass_LiquidSpace += FW_mass_flow_ks_sec
        
        
    def RemoveSteamFlow(self, steamDemand):
        self.runSteamflow = True
        if self.Steam_Mass < 1:
            return
        if self.Steam_Mass < steamDemand:
            steamDemand = self.Steam_Mass
        self.Steam_Mass -= steamDemand
        self.LastSteamDemand = steamDemand
        self.idealGas_PressurePa()

    def shink_swell(self, state):
        if state.region != 4:
            return
        VaporFlow = self.LastSteamDemand * state.Vapor.v
        Sup_Flow = VaporFlow / self.Area
        BRTop = 1.53*math.pow(state.sigma*self.gravity*(state.Liquid.rho-state.Vapor.rho),.25)
        BRBot = math.pow(state.Liquid.rho,0.5)
        BubbleRise = BRTop/BRBot
        RadialDist = 1.2 - 0.2 * math.sqrt(state.Vapor.rho / state.Liquid.rho)
        void_swell = Sup_Flow / (2*BubbleRise + RadialDist*Sup_Flow)
        Swell_Height = self.LiquidHeight / (1-void_swell)
        Swell_Height /= 100
        self.LastSwellHeight = Swell_Height

    def run_Step(self):

        if self.runFW == False or self.runSteamflow == False:
            return
        self.RunFW = False
        self.runSteamflow = False

        self.Enthalpy += self.heatingPower * self.deltaT
        state = iapws.IAPWS97(P=self.PressureMPa, h=self.Enthalpy)

        if state.T != self.TempK:
            self.Enthalpy += ((self.InternalMass * self.carbonSteelCp * (self.TempK - state.T))/self.InternalMass) * self.deltaT
            state = iapws.IAPWS97(P=self.PressureMPa, h=self.Enthalpy)

        self.TempK = state.T
        self.Quality = state.x

        self.Liquid_Mass = self.Total_Mass_LiquidSpace * (1-state.x)
        self.DropletMass = self.Total_Mass_LiquidSpace * state.x

        if (self.DropletMass > 1 and state.region == 4):
            self.Droplet_Travel(state)

        self.WaterVol = self.Liquid_Mass * state.Liquid.v
        if (state.region == 4 and self.DropletMass > 1):
            self.DropletVol = self.DropletMass * state.Vapor.v 
            #self.DropletVol = self.DropletMass * state.Liquid.v #TODO: This is all messed up
            #self.DropletVol = self.DropletMass * state.v #TODO: This is all messed up
        else:
            self.DropletVol = 0
        self.LiquidSpaceVol = self.WaterVol + self.DropletVol
        self.SteamVol = self.Volume - self.LiquidSpaceVol
        self.LiquidHeight = self.LiquidSpaceVol / self.Area

        self.shink_swell(state)
        self.LastLiquidVol = self.LiquidSpaceVol
        self.LiquidHeight += self.LastSwellHeight
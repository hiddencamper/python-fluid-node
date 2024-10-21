import iapws
import iapws._iapws97Constants
import numpy as np
import matplotlib.pyplot as plt
import math

class PID_Controller:
    def __init__(self, setP, setI, setD, target, clampmin, clampmax,inv):
        self.P = setP
        self.I = setI
        self.D = setD
        self.setPoint = target
        self.prev_Int = 0
        self.prev_err = 0
        self.min = clampmin
        self.max = clampmax
        self.inverse = inv

    def run(self, value, dt):
        if self.inverse == False:
            error = self.setPoint - value
        else:
            error = value - self.setPoint
        integral = self.prev_Int + error*dt
        derivative = (error - self.prev_err)/dt
        proportional = error
        self.prev_err = error
        self.prev_Int = integral
        output = self.P * proportional + self.I * integral + self.D * derivative
        
        if output > 100 or output < 0:
            self.prev_Int -= self.I*self.prev_Int

        if output < self.min:
            output = self.min
        if output > self.max:
            output = self.max

        return output/100 #Percentage output


PSI_TO_PASCAL = 6895
FEEDPUMP_MAX_FLOW = 100000 #kg/sec
BYPASS_VLV_FLOW = 60000 #kg/sec
gasConstant = 8.3145
molarMassOfSteam = 18.01528 / 1000 # kg/ mol


Vessel_Area = math.pi * math.pow((251 * 2.54/100)/2,2)
Vessel_Volume = 100000 #m3
Liquid_Volume = 1000 #m3
Temp = 50 + 273.15 #C
PressMPa = 14.7 * PSI_TO_PASCAL / 1e6 #MPa
PressPa = PressMPa * 1e6

state = iapws.IAPWS97(T=Temp, P=PressMPa)
Enthalpy = state.h
Liquid_Mass = Liquid_Volume * state.rho #kg

Steam_Mass = (PressPa*(Vessel_Volume-Liquid_Volume)*molarMassOfSteam) /(Temp * gasConstant)

print("Initial Liquid Mass: ", Liquid_Mass)
print("Initial Steam Mass: ", Steam_Mass)

steps = 1000
stepsize=1

WtrMassArr = np.zeros(steps)
StmMassArr = np.zeros(steps)
EnthalpyArr = np.zeros(steps)
TempArr = np.zeros(steps)
PressArr = np.zeros(steps)
WaterHeightArr = np.zeros(steps)
QualityArr = np.zeros(steps)
FeedwaterFlowArr = np.zeros(steps)
HeatArr = np.zeros(steps)
SteamFlowArr = np.zeros(steps)

heat_add = 0
heat_limit = 100
PressPSIG = 0

FW_PID = PID_Controller(4, .5, .25,30,0,100,False)

SB_PID = PID_Controller(1,.5,.2,1000,0,100, True)

for i in range (steps):
    try:
        Enthalpy += heat_add
        print(i)

        state = iapws.IAPWS97(P=PressMPa, h=Enthalpy)

#        if i == 471:
#            print ("Stop")

        TempArr[i] = state.T-273.15
        Temp = TempArr[i]
        WtrMassArr[i] = Liquid_Mass * (1-state.x)

        StmMassArr[i] = Steam_Mass + Liquid_Mass*state.x
        EnthalpyArr[i] = Enthalpy
        QualityArr[i] = state.x
        SteamFlowArr[i] = 0
        SteamFlowArr[i] = BYPASS_VLV_FLOW * SB_PID.run(PressPSIG,1)

        if i > 300 and i < 700:
            SteamFlowArr[i] = StmMassArr[i] * .1
            if SteamFlowArr[i] > BYPASS_VLV_FLOW:
                SteamFlowArr[i] = BYPASS_VLV_FLOW

        if SteamFlowArr[i] > 0 and state.region == 4:
            Enthalpy -= (SteamFlowArr[i] * state.Vapor.h) / (StmMassArr[i] + WtrMassArr[i])

        StmMassArr[i] -= SteamFlowArr[i]
        Liquid_Mass -= SteamFlowArr[i]
        Liquid_Volume = WtrMassArr[i] * state.Liquid.v
        WaterHeightArr[i] = Liquid_Volume / Vessel_Area
        
        Steam_Volume = Vessel_Volume - Liquid_Volume
 
        PressPa = (StmMassArr[i]/molarMassOfSteam) * gasConstant * state.T / Steam_Volume # Pa
        PressMPa = PressPa/1e6
        PressPSIG = PressPa/ PSI_TO_PASCAL

        new_state = iapws.IAPWS97(P=PressMPa, h=Enthalpy)


        PressArr[i] = PressPSIG
        if PressPSIG >= 1100:
            heat_add -=5
            if heat_add < 0:
                heat_add = 0
        else:
            heat_add += 5
            if heat_add > heat_limit:
                heat_add = heat_limit

        if i > 500 and i < 600:
           heat_add = 0
        HeatArr[i] = heat_add



        FeedWaterKg = FEEDPUMP_MAX_FLOW * FW_PID.run(WaterHeightArr[i],1)

        FeedwaterFlowArr[i] = FeedWaterKg
        Liquid_Mass += FeedWaterKg
        Enthalpy = (WtrMassArr[i] * EnthalpyArr[i] + FeedWaterKg*200) / (WtrMassArr[i] + FeedWaterKg)
        EnthalpyArr[i] = Enthalpy
        WtrMassArr[i] += FeedWaterKg
        Liquid_Volume= WtrMassArr[i] * state.Liquid.v        


    except:
        break


time_array = np.arange(0,steps, 1)
plt.figure(figsize=(16, 8))

plt.subplot(5, 1, 1)
plt.plot(time_array, SteamFlowArr)
plt.xlabel("Time (s)")
plt.ylabel("Steam Flow - SRV (kg/sec)")
plt.title("Steam Flow Over Time")

plt.subplot(5, 1, 2)
plt.plot(time_array, PressArr)
plt.xlabel("Time (s)")
plt.ylabel("Pressure (PSIG)")
plt.title("Pressure Over Time")

plt.subplot(5, 1, 3)
plt.plot(time_array, EnthalpyArr)
plt.xlabel("Time (s)")
plt.ylabel("Sp Enthalpy (kj/kg)")
plt.title("Sp Enthalpy Over Time")

plt.subplot(5, 1, 4)
plt.plot(time_array, WaterHeightArr)
plt.xlabel("Time (s)")
plt.ylabel("Water Height")
plt.title("Water Height Over Time")

plt.subplot(5, 1, 5)
plt.plot(time_array, FeedwaterFlowArr)
plt.xlabel("Time (s)")
plt.ylabel("Feedwater (kg/sec)")
plt.title("Feedwater Flow Over Time")
    
plt.tight_layout()
plt.show()
    


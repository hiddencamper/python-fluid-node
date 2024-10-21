import iapws
import iapws._iapws97Constants
import numpy as np
import matplotlib.pyplot as plt
import math
import PID
import boiler_func

PSI_TO_PASCAL = 6895
FEEDPUMP_MAX_FLOW = 100000 #kg/sec
BYPASS_VLV_FLOW = 60000 #kg/sec
gasConstant = 8.3145
molarMassOfSteam = 18.01528 / 1000 # kg/ mol

Vessel_Diameter = 218 * 2.54/100 # inches converted to cm converted to m
Vessel_Radius = Vessel_Diameter/2
Vessel_Area = math.pi * Vessel_Radius * Vessel_Radius
Vessel_Volume = 100000 #m3
Vessel_Height = Vessel_Volume / Vessel_Area
Liquid_Volume = 1000 #m3
Temp = 50 + 273.15 #C
PressMPa = 14.7 * PSI_TO_PASCAL / 1e6 #MPa
PressPa = PressMPa * 1e6
Internal_Mass = 10000 #kg

steps = 500
stepsize= 1

Boiler = boiler_func.Cylindrical_Boiler(Vessel_Radius, Vessel_Height, Liquid_Volume, Temp, PressMPa, Internal_Mass,stepsize)

print("Initial Liquid Mass: ", Boiler.Liquid_Mass)
print("Initial Steam Mass: ", Boiler.Steam_Mass)



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
DropletMassArr = np.zeros(steps)
DropletFlowArr = np.zeros(steps)
DropletVelocityArr = np.zeros(steps)

heat_add = 0
heat_limit = 50
Boiler.Set_Heating_Power(heat_add)

FW_PID = PID.PID_Controller(4, 1, .25,40,0,100,False)

SB_PID = PID.PID_Controller(1,.5,.2,1000,0,100, True)

for i in range (steps):
    try:
        print(i)

        if i == 31:
            print("stop")

        Boiler.Set_Heating_Power(heat_add)

        Boiler.run_Step()

        FeedWaterKg = FEEDPUMP_MAX_FLOW * FW_PID.run(Boiler.LiquidHeight,stepsize)

        Boiler.injectFW(FeedWaterKg, 200)

        SteamFlowKg = BYPASS_VLV_FLOW * SB_PID.run(Boiler.PressPSIG,stepsize)

        Boiler.RemoveSteamFlow(SteamFlowKg)

        TempArr[i] = Boiler.TempK
        WtrMassArr[i] = Boiler.Liquid_Mass
        DropletMassArr[i] = Boiler.DropletMass
        StmMassArr[i] = Boiler.Steam_Mass
        EnthalpyArr[i] = Boiler.Enthalpy
        QualityArr[i] = Boiler.Quality
        SteamFlowArr[i] = SteamFlowKg
        WaterHeightArr[i] = Boiler.LiquidHeight
        PressArr[i] = Boiler.PressPSIG
        HeatArr[i] = Boiler.heatingPower
        FeedwaterFlowArr[i] = FeedWaterKg
        DropletFlowArr[i] = Boiler.LastDropletFlow
        DropletVelocityArr[i] = Boiler.LastDropletVel

        if Boiler.PressPSIG >= 1100:
            heat_add -=5
            if heat_add < 0:
                heat_add = 0
        else:
            heat_add += 5
            if heat_add > heat_limit:
                heat_add = heat_limit 


    except:
        break


time_array = np.arange(0,steps, 1)
plt.figure(figsize=(16, 8))

plt.subplot(5, 1, 1)
plt.plot(time_array, SteamFlowArr)
plt.xlabel("Time (s)")
plt.ylabel("Steam Flow (kg/s)")
plt.title("Steam Flow Over Time")

plt.subplot(5, 1, 2)
plt.plot(time_array, PressArr)
plt.xlabel("Time (s)")
plt.ylabel("Pressure (PSIG)")
plt.title("Pressure Over Time")

plt.subplot(5, 1, 3)
plt.plot(time_array, QualityArr)
plt.xlabel("Time (s)")
plt.ylabel("Quality")
plt.title("Quality Over Time")

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
    


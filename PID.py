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
        self.ThreeElement = False        
        self.ThreeElementKp = 0

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

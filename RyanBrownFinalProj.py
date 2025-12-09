import math, matplotlib.pyplot as plt, scipy, numpy as np

class Motor():
    def __init__(self, motorChoice):
        self.name = motorChoice.upper()

        if self.name == "B6": #https://www.thrustcurve.org/simfiles/5f4294d20002e90000000414/
            self.thrustCurve = [
                (0.036, 1.364),
                (0.064, 2.727),
                (0.082, 4.215),
                (0.111, 6.694),
                (0.146, 9.545),
                (0.172, 11.901),
                (0.181, 12.149),
                (0.191, 11.901),
                (0.211, 9.174),
                (0.239, 7.314),
                (0.264, 6.074),
                (0.275, 5.95),
                (0.333, 5.207),
                (0.394, 4.835),
                (0.445, 4.835),
                (0.556, 4.339),
                (0.667, 4.587),
                (0.723, 4.339),
                (0.793, 4.091),
                (0.812, 2.603),
                (0.833, 1.24),
                (0.857, 0)
            ]
            self.propMass = 0.006
            self.propBurnTime = 0.8
        elif self.name == "E30":
            self.thrustCurve = [
                (0.01, 49),
                (0.02, 49),
                (0.05, 46),
                (0.10, 44),
                (0.20, 43),
                (0.25, 42),
                (0.30, 41),
                (0.35, 40),
                (0.40, 39),
                (0.45, 38),
                (0.50, 37),
                (0.55, 35),
                (0.60, 33),
                (0.65, 32),
                (0.70, 31),
                (0.75, 30),
                (0.80, 27),
                (0.85, 25),
                (0.90, 20),
                (0.91, 19),
                (0.93, 12),
                (0.95, 5),
                (0.97, 1),
                (1.00, 0)
            ]
            self.propMass = 0.018
            self.propBurnTime = 1.0
        elif self.name == "G80":
            self.thrustCurve = [
                (0.013, 89.054),
                (0.018, 101.584),
                (0.029, 105.388),
                (0.047, 102.927),
                (0.104, 100.018),
                (0.19, 102.255),
                (0.268, 104.94),
                (0.306, 104.269),
                (0.352, 105.835),
                (0.41, 104.94),
                (0.432, 106.73),
                (0.512, 105.612),
                (0.563, 106.954),
                (0.591, 104.493),
                (0.62, 106.283),
                (0.715, 103.15),
                (0.76, 103.374),
                (0.806, 101.36),
                (0.848, 102.032),
                (1.2, 88.383),
                (1.243, 77.866),
                (1.269, 57.952),
                (1.302, 44.079),
                (1.329, 34.011),
                (1.367, 26.403),
                (1.398, 22.823),
                (1.433, 20.585),
                (1.48, 21.704),
                (1.508, 20.809),
                (1.584, 13.201),
                (1.604, 11.635),
                (1.652, 11.188),
                (1.683, 5.37),
                (1.701, 0)
            ]
            self.propMass = 0.063
            self.propBurnTime = 1.7

        elif self.name == "J250":
            self.thrustCurve = [
                (0.005, 149.863),
                (0.016, 124.444),
                (0.048, 224.0),
                (0.051, 290.723),
                (0.077, 309.787),
                (0.221, 324.085),
                (0.269, 316.141),
                (0.497, 335.735),
                (0.997, 324.085),
                (1.497, 288.605),
                (1.713, 266.893),
                (2.003, 210.231),
                (2.162, 176.87),
                (2.487, 119.678),
                (2.516, 125.503),
                (2.718, 91.612),
                (2.766, 45.541),
                (2.814, 18.005),
                (2.867, 7.943),
                (2.902, 0.0)
            ]
            self.propMass = 0.487
            self.propBurnTime = 2.9

        else:
            raise ValueError("Unknown Motor Type")
        
        self.timePoints = np.array([p[0] for p in self.thrustCurve]) #creates an array for the thrust curve, one for thrust, and another for its associated time during burnout
        self.thrustPoints = np.array([p[1] for p in self.thrustCurve])
    
    def getThrust(self, t):
        #returns thrust in Newtons, given time t
        return float(np.interp(t, self.timePoints, self.thrustPoints, left=0.0, right = 0.0)) #uses NumPy interpolation for any time t, to find an accurate thrust estimate

class Rocket():
    def __init__(self, bodyMass, motorChoice, Cd, area, wind, chuteHeight, chuteArea, dt):
        self.mass = bodyMass #mass of rocket structure, without fuel
        self.motor = Motor(motorChoice) #motor choice + creates motor object
        self.Cd = Cd #coefficient of drag
        self.area = area #cross-sectional area of rocket
        self.wind = wind #horizontal wind velocity (positive is to right)
        self.chuteHeight = chuteHeight #height which chute is deployed to decelerate rocket
        self.chuteArea = chuteArea #size of chute to determine deceleration
        self.timeStep = dt #time step for simulation, smaller = more precise but slower simulation
        self.x = [0.0] #position of rocket
        self.v = [0.0] #velocity of rocket
        self.a = [0.0] #acceleration of rocket
        self.landed = False
        self.apogee = 0.0 #initial value for apogee
        self.chuteDeploy = False #boolean expression to model whether the chute has been deployed
        self.chuteDeployTime = 0.0
        print("\n \n Rocket created! \n")

    def getThrust(self, t):
        return self.motor.getThrust(t) #uses Motor method to find thrust at time t
    
    def getDrag(self, alt, vel, t,):
        rho = atm(alt)[2]
        if self.chuteDeploy == False:
            drag = .5*rho*vel*vel*self.Cd*self.area #calculates the drag force on the rocket at a given velocity and altitude
        else:
            chuteCd = 1.5*(1-np.exp(-0.4*(t-self.chuteDeployTime)))
            #print(f"Time: {t}, chuteCD: {chuteCd}")
            drag = .5*rho*vel*vel*self.Cd*self.area + .5*rho*vel*vel*chuteCd*self.chuteArea #calculates drag, including the drag from the chute, using a default Cd for chutes of 1.5 across all sims
        drag *= -np.sign(vel) #flips the direction of drag, so it is always in the opposite direction of velocity
        return drag


def atm(alt): 
    """
    Function to calculate atmosphere conditions, including temperature, pressure, 
    and density based on ISA Equations
    """
    # T --> Kelvin, rho --> kg/m^3, p --> Pa

    h = max(alt, 0) #troubleshoots problem at beginning of sim where h sometimes becomes negative

    if h < 11000:  # Troposphere
        T = 288.15 - 0.0065 * h
        p = 101325 * (T / 288.15) ** 5.2559

    elif 11000 <= h < 25000:  # Lower Stratosphere
        T = 216.65
        p = 22632.06 * math.exp(-9.80665 * (h - 11000) / (287.05 * T))

    else: #h >= 25000:  # Mid Stratosphere
        T = 216.65 + 0.00299 * (h - 25000)
        p = 2488.0 * (T / 216.65) ** (-9.80665 / (287.05 * 0.00299))

    rho = p / (287.05 * T)
    return T, p, rho

def grav(h): 
    "Calculates gravitational acceleration at altitude h"
    g = 9.81*(1+(h/(6371*10**3)))**(-2)

    return g

def export_data(file, rocket, time): #function to export the simulation data to a .txt file
    with open(file, 'w') as f:
        f.write("#Rocket Parameters \n")
        f.write(f"Rocket Mass: {rocket.mass} kg\n")
        f.write(f"Coefficient of Drag: {rocket.Cd}\n")
        f.write(f"Rocket Face Area: {rocket.area} m^2\n")
        f.write(f"Motor Choice: {rocket.motor.name}\n")
        f.write(f"Propellant Mass: {rocket.motor.propMass} kg\n")
        f.write(f"Motor Burn Time: {rocket.motor.propBurnTime} s\n")
        f.write(f"Chute Deployment Height: {rocket.chuteHeight} m\n")
        f.write(f"Time Step (dt): {rocket.timeStep} s\n\n")

        f.write('# Flight Data (time, altitude, velocity, acceleration)\n')
        for i in range(len(time)):
            f.write(f"{time[i]}, {rocket.x[i]}, {rocket.v[i]}, {rocket.a[i]}\n")

    print(f"\nSaved simulation details to {file}\n")

def setup():
    while True:
        choice = input("Would you like to choose a preset rocket, or use custom parameters? (preset or custom) ")

        if choice.lower() == "custom":
            bodyMass = float(input("What is the mass of the rocket, excluding the mass of the fuel? (kg) "))
            motorChoice = input("What is your choice of motor? (B6, E30, G80, J250) ")
            coefDrag = float(input("What is the coefficient of drag for your rocket? "))
            area = float(input("What is the cross-sectional area of your rocket? (m^2) "))
            wind = float(input("What is the velocity of the wind? (Positive is to right) (m/s) "))
            chuteHeight = float(input("At what altitude will the chute deploy? (m) "))
            chuteArea = float(input("What is the cross-sectional area of your chute? (m^2) "))
            timeStep = float(input("What is the time-step for this simulation? (Smaller = more precision) "))
            return Rocket(bodyMass, motorChoice, coefDrag, area, wind, chuteHeight, chuteArea, timeStep)
        
        elif choice.lower() == "preset":
            print("""
                  Choice A: Alpha Trainer (m=.25 kg, B6 Motor, Cd = .55, A = .0025 m^2)

                  Choice B: Sky Piercer (m = 1.1 kg, E30 motor, Cd = .5, Area = .0038 m^2)

                  Choice C: High Flyer (m = 1.8 kg, G80 motor, Cd = .45, Area = .005 m^2)

                  Choice D: Chonk Lifter (m = 4 kg, J250 motor, Cd = .6, Area = .0075 m^2)

                  """)
            presetChoice = input("Which preset would you like to choose? (A, B, C, or D) ")

            if presetChoice.upper() == "A":
                return Rocket(.225, "B6", .55, .0025, 0, 10, .9, .001)
            elif presetChoice.upper() == "B":
                return Rocket(1.1, "E30", .45, .0038, 0, 25, 1.55, .001)
            elif presetChoice.upper() == "C":
                return Rocket(1.8, "G80", .35, .005, 0, 75, 2.25, .001)
            else:
                return Rocket(4, "J250", .55, .0075, 0, 150, 3.15, .001)
            
        else:
            print("Invalid input, try again. \n \n")



def main():
    rocket = setup()

    dt = rocket.timeStep

    t = [0.0]
    force = [rocket.getThrust(t[0]) + rocket.getDrag(0, 0, 0) - (rocket.mass+ rocket.motor.propMass)*grav(rocket.x[0])]
    rocket.a[0] = force[0]/(rocket.mass + rocket.motor.propMass)
    i = 1

    plt.ion()
    fig, ax = plt.subplots()
    line, = ax.plot([], [])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Altitude (m)")
    ax.set_title("Live Altitude")
    plt.show()

    while rocket.landed == False: #runs the loop until the rocket has landed (returned to x = 0)
        
        t.append(dt*i) #creates time list, with each element being the last element + timestep

        if t[i] < rocket.motor.propBurnTime:
            mass = rocket.mass + rocket.motor.propMass*(1-t[i]/rocket.motor.propBurnTime)
        else:
            mass = rocket.mass

        force.append(rocket.getThrust(t[i]) + rocket.getDrag(rocket.x[i-1], rocket.v[i-1], t[i]) - mass*grav(rocket.x[i-1])) #calulates the force on the rocket at time t

        rocket.a.append(force[i]/mass) #calculates the acceleration of the rocket at this increment, using Newton's second law
        rocket.v.append(rocket.v[i-1] + rocket.a[i]*dt)
        rocket.x.append(rocket.x[i-1] + rocket.v[i-1]*dt + .5*rocket.a[i-1]*dt*dt)

        if rocket.x[i] > rocket.apogee:
            rocket.apogee = rocket.x[i]
            apogee_index = i
            apogee_time = dt*i


        if rocket.apogee > rocket.chuteHeight and rocket.x[i] <= rocket.chuteHeight and rocket.chuteDeploy == False:
            rocket.chuteDeploy = True
            rocket.chuteDeployTime = t[i]

        if rocket.x[i] <= 0 and i > 500: #i > 500 to ensure the simulation runs for at least 1/2 second, negating any interpolation errors
            rocket.x[i] = 0
            rocket.landed = True
        if i % 100 == 0: #update the plot every 100 iterations
            line.set_xdata(t)
            line.set_ydata(rocket.x)
            ax.relim()
            ax.autoscale_view()
            plt.pause(dt)
            #print("ALT:", rocket.x[i])

        i += 1
        
    peak_vel = max(rocket.v)
    peak_vel_index = rocket.v.index(peak_vel)
    peak_vel_time = dt*peak_vel_index

    plt.close(fig) #closes the interactive plots
    plt.ioff() #disables itneractive mode

    plt.subplot(2,1,1)
    plt.plot(t, rocket.x)
    plt.scatter(apogee_time, rocket.apogee, color='red') # type: ignore
    plt.text(apogee_time, rocket.apogee, f"Apogee: {rocket.apogee:.1f} m @ {apogee_time:.2f} s") # type: ignore
    plt.title("Altitude vs Time")
    plt.xlabel('Time (s)')
    plt.ylabel('Altitude (m)')

    plt.subplot(2,1,2)
    plt.plot(t, rocket.v)
    plt.scatter(peak_vel_time, peak_vel, color = 'red')
    plt.text(peak_vel_time, peak_vel, f"Peak v: {peak_vel:.1f} m/s @ {peak_vel_time:.2f} s")
    plt.title('Velocity vs Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity (m/s)')
    plt.tight_layout()
    plt.show()

    fileChoice = input("Would you like to create a .txt file containing all the simulation data? (Yes or No) ")
    if fileChoice.lower() == "yes":
        fileName = input("What would you like to name the file? ") + ".txt"
        export_data(fileName, rocket, t)
    


main()


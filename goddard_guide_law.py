import numpy as np
import matplotlib.pyplot as plt


class Rocket:
    def __init__(self):
        self.mass = 17500.              # kg
        self.thrust_max = 20.*17500        # N
        self.thrust = 0.            # N
        self.mass_min = 10000.          # kg
        self.height = 0.            # m
        self.gravity = 9.81         # m/s^2
        self.velocity = 0.
        self.specific_impulse = 330.  # s
        self.drag = np.zeros([4, 1])              # N
        self.alpha = 1./self.gravity/self.specific_impulse
        self.c = 1./self.alpha
        self.time = 0.
        self.dt = 1.              # time step, 60 Hz
        self.is_sing = False
        self.sing_tol = 1e-6
        self.state = np.zeros([4, 1])

    def rocket_state(self, t, s):
        # rocket state equations
        # s = [h, v, m]
        # integrate height over timestep..
        snew = np.array([0., 0., 0., 0.])

        # get thrust
        thrust = self.get_thrust(s)
        snew[0] = s[1]
        d = self.drag[0]
        snew[1] = thrust/s[2] - d/s[2] - self.gravity
        snew[2] = -self.alpha*thrust
        snew[3] = thrust
        return snew

    def rk4(self, t, s):
        # RK4 Algorithm:
        # y_(n+1) = y_n + 1/6*(k_1 + 2*k_2 + 2*k_3 + k_4)*h
        # t_(n+1) = t_n + h
        # k_1 = f(t_n, y_n)
        # k_2 = f(t_n + h/2, y_n + h*k_1/2)
        # k_3 = f(t_n + h/2, y_n + h*k_2/2)
        # k_4 = f(t_n + h, y_n + h*k_3)

        k1 = self.rocket_state(t, s)
        k2 = self.rocket_state(t + self.dt*0.5, s + self.dt*k1*0.5)
        k3 = self.rocket_state(t + self.dt*0.5, s + self.dt*k2*0.5)
        k4 = self.rocket_state(t + self.dt, s + self.dt*k3)
        self.time = t + self.dt
        return s + 1./6.*(k1 + 2.*k2 + 2.*k3 + k4)*self.dt

    # is the singularity thrust to be calculated?
    def check_singularity(self):
        self.is_sing = self.sing_tol >= \
                        self.drag[0] + self.state[2]*self.gravity - \
                        self.alpha*self.drag[0]*self.state[1] - \
                        self.state[1]*self.drag[1]

    def get_thrust(self, s):
        #self.check_singularity(s)
        if self.is_sing:
            # update drag model.
            #drag = self.calculate_drag(s)
            d = self.drag[0]
            dv = self.drag[1]
            ddv = self.drag[2]
            ddh = self.drag[3]
            ddvdh = 0.  # update from drag terms with updated model
            p1 = s[2]/(d + 2.*self.c*dv+self.c**2*ddv)
            p2 = -self.gravity*(d+self.c*dv) #+self.c*(self.c-s[1])*ddh-s[1]*self.c**2*ddvdh
            #thrust = d + s[2]*self.gravity + p1*p2
            thrust = d + s[2]*self.gravity + s[2]/(d+2*self.c*dv+self.c**2*ddv)*-self.gravity*(d+self.c*dv)
            thrust = float(thrust)  # remove warning 
        else:
            thrust = self.thrust_max
        if thrust > self.thrust_max:
            raise Exception("THRUST CANNOT BE BIGGER THAN THRUST_MAX!!!")
        return thrust

    def calculate_drag(self):
        # simple drag model that should be updated for complexity--this mirrors the model used in ECE6570
        drag =  np.array([0., 0., 0., 0.])
        drag[0] = 0.5*self.state[1]**2
        drag[1] = self.state[1]
        drag[2] = 1.
        drag[3] = self.state[0]
        return drag

    def set_state(self):
        self.state[0] = self.height
        self.state[1] = self.velocity
        self.state[2] = self.mass
        self.state[3] = self.thrust

    def set_state_params(self):
        self.height = self.state[0]
        self.velocity = self.state[1]
        self.mass = self.state[2]
        self.thrust = self.state[3]

    # schedule updates to state, drag, etc.
    def update(self):
        # update drag and singularity calculations first
        self.drag = self.calculate_drag()
        # check singularity
        self.check_singularity()
        self.height, self.velocity, self.mass, dummy = self.rk4(self.time, np.array([self.height, self.velocity, self.mass, self.thrust]))
        self.thrust = self.get_thrust(self.state) # uncumbered by rk4 algorithm
        # update state params
        self.set_state()



def main():
    rocket = Rocket()
    rocket.set_state() # initial parameters set to state
    tsav = []
    thsav = []
    hsav = []
    sing = []
    while rocket.mass >= rocket.mass_min:
        # print(rocket.mass)
        rocket.update()
        # print(rocket.mass)
        tsav.append(rocket.time)
        thsav.append(rocket.thrust)
        hsav.append(rocket.height)
        sing.append(rocket.is_sing)
    np.array(tsav)
    np.array(thsav)
    np.array(hsav)
    np.array(sing)
    plt.figure("Thrust v time")
    plt.plot(tsav, thsav)
    plt.show()
    plt.figure("height v time")
    plt.plot(tsav, hsav)
    plt.show()

if __name__ == '__main__':
    main()


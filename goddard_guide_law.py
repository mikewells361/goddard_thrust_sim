import numpy as np


class Rocket:
    def __init__(self):
        self.mass = 0.              # kg
        self.thrust_max = 0.        # N
        self.thrust = 0.            # N
        self.mass_init = 0.         # kg
        self.mass_min = 0.          # kg
        self.height = 0.            # m
        self.height_init = 0.       # m
        self.gravity = 9.81         # m/s^2
        self.velocity = 0.
        self.velocity_init = 0.     # m/s
        self.specific_impulse = 100.  # s
        self.drag = np.zeros([4, 1])              # N
        self.alpha = 1./self.gravity/self.specific_impulse
        self.c = 1./self.alpha
        self.time = 0.
        self.is_sing = False
        self.sing_tol = 1e-6


    def rocket_state(self):
        pass

    def rk4(self):
        pass

    # is the singularity thrust to be calculated?
    def check_singularity(self):
        self.calculate_drag()
        self.sing_tol = self.sing_tol >= \
                        self.drag[0] + self.mass*self.gravity - \
                        self.alpha*self.drag[0]*self.velocity - \
                        self.velocity*self.drag[1]

    def get_thrust(self):
        self.check_singularity()
        if self.is_sing:
            # update drag model. also drag is updated in check_singularity so no need to update it again here
            d = self.drag[0]
            dv = self.drag[1]
            ddv = self.drag[2]
            ddh = self.drag[3]
            ddvdh = 0.  # update from drag terms with updated model

            p1 = self.mass/(d + 2.*self.c*dv+self.c**2*ddv)
            p2 = -self.gravity*(d+self.c*dv)+self.c*(self.c-self.velocity)*ddh-self.velocity*self.c**2*ddvdh
            self.thrust = d + self.mass*self.gravity + p1*p2
        else:
            self.thrust = self.thrust_max

    def calculate_drag(self):
        # simple drag model that should be updated for complexity--this mirrors the model used in ECE6570
        self.drag[0] = 0.5*self.velocity**2
        self.drag[1] = self.velocity
        self.drag[2] = 1.
        self.drag[3] = self.height
        pass

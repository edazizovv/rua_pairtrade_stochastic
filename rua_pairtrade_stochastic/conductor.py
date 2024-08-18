import numpy


class Governor:
    def __init__(self, T, gamma, X_k, X_theta, X_eta, B_sigma, BX_ZW_rho):
        self.T = T
        self.gamma = gamma
        self.k = X_k
        self.theta = X_theta
        self.eta = X_eta
        self.sigma = B_sigma
        self.rho = BX_ZW_rho
        self.x = numpy.nan
        self.t = numpy.nan

    def step(self, x, t):
        self.x = x
        self.t = t
        result = self.compute_h()
        return result

    def compute_h(self):
        result = self._h_1() * (self._h_2() + self._h_3() + self._h_4())
        return result

    def compute_alpha(self):
        result = self._a_1() * self._a_2()
        return result

    def compute_beta(self):
        result = self._b_1() * (self._b_2() + self._b_3())
        return result

    def _h_1(self):
        result = 1 / (1 - self.gamma)
        return result

    def _h_2(self):
        result = self.compute_beta() + 2 * self.x * self.compute_alpha()
        return result

    def _h_3(self):
        result = -1 * self.k * (self.x - self.theta) / (self.eta ** 2)
        return result

    def _h_4(self):
        result = (self.rho * self.sigma / self.eta) + 0.5
        return result

    def _a_1(self):
        result = (self.k * self._u_1()) / (2 * self.eta ** 2)
        return result

    def _a_2(self):
        result = 1 + (2 * numpy.sqrt((1 - self.gamma))) / (self._u_1() - self._u_2() * self._u_3())
        return result

    def _b_1(self):
        result = 1 / (2 * self.eta ** 2 * (self._u_1() + self._u_2() * self._u_3()))
        return result

    def _b_2(self):
        result = self.gamma * numpy.sqrt((1 - self.gamma)) * self._u_4() * (1 - self._u_3()) ** 2
        return result

    def _b_3(self):
        result = -1 * self.gamma * (self._u_4() + 2 * self.k * self.theta) * (1 - self._u_3())
        return result

    def _u_1(self):
        result = 1 - numpy.sqrt((1 - self.gamma))
        return result

    def _u_2(self):
        result = 1 + numpy.sqrt((1 + self.gamma))
        return result

    def _u_3(self):
        result = numpy.exp(((2 * self.k * (self.T - self.t)) / numpy.sqrt((1 - self.gamma))))
        return result

    def _u_4(self):
        result = self.eta ** 2 + 2 * self.rho * self.sigma * self.eta
        return result


def govern():
    return None

def governor(gamma):
    # takes parameters == returns weights

    pre_multiplier = 1 / (1 - gamma)
    beta_first
    return None
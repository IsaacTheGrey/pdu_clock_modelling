from scipy.integrate import odeint
from .model import goodwin_model_simplified

def integrate_model(y0, t, parameters):
    """
    Integrates the ODEs numerically.
    """
    return odeint(goodwin_model_simplified, y0, t, args=(parameters,))

from mewpy.simulation import get_container, get_simulator


def load_sbml_container(filename, flavor='reframed'):
    if flavor == 'reframed':
        from reframed.io.sbml import load_cbmodel
        model = load_cbmodel(filename, flavor='cobra')
    elif flavor == 'cobra':
        from cobra.io import read_sbml_model
        model = read_sbml_model(filename)
    else:
        raise ValueError(f"{flavor} is not a recognized flavor")
    container = get_container(model)
    return container


def load_sbml_simulator(filename, flavor='reframed', envcond=None):
    if flavor == 'reframed':
        from reframed.io.sbml import load_cbmodel
        model = load_cbmodel(filename, flavor='cobra')
    elif flavor == 'cobra':
        from cobra.io import read_sbml_model
        model = read_sbml_model(filename)
    else:
        raise ValueError(f"{flavor} is not a recognized flavor")
    simul = get_simulator(model, envcond=envcond)
    return simul


def load_gecko_simulator(filename, flavor='reframed'):
    if flavor == 'reframed':
        from mewpy.model.gecko import GeckoModel
        model = GeckoModel(filename)
    elif flavor == 'cobra':
        from geckopy.gecko import GeckoModel
        model = GeckoModel(filename)
    else:
        raise ValueError(f"{flavor} is not a recognized flavor")
    simul = get_simulator(model)
    return simul

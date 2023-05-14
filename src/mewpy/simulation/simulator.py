
# Copyright (C) 2019- Centre of Biological Engineering,
#     University of Minho, Portugal

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
##############################################################################
Author: Vítor Pereira
##############################################################################
"""
from .simulation import Simulator
from mewpy.util.constants import ModelConstants

# Model specific simulators mapping:
# Entries take the form:  full_model_class_path -> (simulator_path, simulator_class_name)
# TODO: use qualified names
map_model_simulator = {
    'geckopy.gecko.GeckoModel': ('mewpy.simulation.cobra', 'GeckoSimulation'),
    'mewpy.model.gecko.GeckoModel': ('mewpy.simulation.reframed', 'GeckoSimulation'),
    'mewpy.model.smoment.SMomentModel': ('mewpy.simulation.reframed', 'GeckoSimulation'),
    'mewpy.germ.models.model.Model': ('mewpy.simulation.germ', 'Simulation'),
    'mewpy.germ.models.metabolic.MetabolicModel': ('mewpy.simulation.germ', 'Simulation'),
    'mewpy.germ.models.regulatory.RegulatoryModel': ('mewpy.simulation.germ', 'Simulation'),
    'mewpy.germ.models.model.MetabolicRegulatoryModel': ('mewpy.simulation.germ', 'Simulation'),
    'mewpy.germ.models.model.RegulatoryMetabolicModel': ('mewpy.simulation.germ', 'Simulation'),
}


def get_simulator(model, envcond=None, constraints=None, reference=None, reset_solver=ModelConstants.RESET_SOLVER):
    """
    Returns a simulator instance for the model.
    The simulator instance is model dependent.
    Besides able to be used on its own, this function is invoked by EA optimization problems
    and by evaluation function instances to perform phenotyoe evaluations of candidate solutions.


    :param model: the model
    :param dict envcond: A dictionary of environmental conditions.
    :param dict constraints: A dictionary of additional persistent constraints.
    :param dict reference: A dictionary of the wild type flux values
    :param bool reset_solver: Whether to reset the solver before each simulation
    :returns: An instance of Simulator
    """
    # already is a Simulator instance
    if isinstance(model, Simulator):
        return model

    instance = None
    name = f"{model.__class__.__module__}.{model.__class__.__name__}"
    if name in map_model_simulator:
        module_name, class_name = map_model_simulator[name]
        module = __import__(module_name, fromlist=[None])
        class_ = getattr(module, class_name)
        try:
            model.solver.configuration.timeout = ModelConstants.SOLVER_TIMEOUT
        except:
            pass
        try:
            model.solver.problem.params.OutputFlag = 0
        except Exception as e:
            pass
        instance = class_(model, envcond=envcond,
                          constraints=constraints, reference=reference, reset_solver=reset_solver)
    elif "etfl" in name:
        try:
            from .cobra import Simulation
            from etfl.optim.config import standard_solver_config
            standard_solver_config(model, verbose=False)
            model.solver.configuration.timeout = max(7200, ModelConstants.SOLVER_TIMEOUT)
            instance = Simulation(
                model, envcond=envcond, constraints=constraints, reference=reference, reset_solver=reset_solver)
            instance._MAX_STR = 'max'
            instance._MIN_STR = 'min'
        except Exception:
            raise RuntimeError("Could not create simulator for the ETFL model")
    else:
        try:
            from cobra.core.model import Model
            if isinstance(model, Model):
                from .cobra import Simulation
                model.solver.configuration.timeout = ModelConstants.SOLVER_TIMEOUT
                instance = Simulation(
                    model, envcond=envcond, constraints=constraints, reference=reference, reset_solver=reset_solver)
        except ImportError:
            pass
        if not instance:
            try:
                from reframed.core.cbmodel import CBModel
                if isinstance(model, CBModel):
                    from .reframed import Simulation
                    instance = Simulation(
                        model, envcond=envcond, constraints=constraints, reference=reference, reset_solver=reset_solver)
            except ImportError:
                pass

    if not instance:
        raise ValueError(f"The model <{name}> has no defined simulator.")
    return instance


def get_container(model):
    """
    Returns a container for a given model sharing a common interface.
    A container does not perform any task, it only serves as a basic interface with a phenotype simulator.

    :param model: An instance of a metabolic model.
    :returns: A container.

    """

    from mewpy.germ.models import Model, RegulatoryModel, MetabolicModel

    if isinstance(model, (Model, MetabolicModel, RegulatoryModel)):

        from .germ import GERMModel

        return GERMModel(model)

    try:
        from reframed.core.cbmodel import CBModel
        if isinstance(model, CBModel):
            from mewpy.simulation.reframed import CBModelContainer
            return CBModelContainer(model)
    except Exception:
        pass

    try:
        from cobra.core.model import Model
        if isinstance(model, Model):
            from mewpy.simulation.cobra import CobraModelContainer
            return CobraModelContainer(model)
    except Exception:
        pass

    raise ValueError(f"Unrecognized model class: {model.__class__.name}")
